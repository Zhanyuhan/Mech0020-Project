//> using scala "3.3"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.92.0"

import scalismo.geometry._
import scalismo.mesh._
import scalismo.io._
import scalismo.utils.Random.FixedSeed.randBasis
import java.io.{File, PrintWriter}
import scala.util.Random

object AnatomicalAlignedToCSV extends App {

  scalismo.initialize()
  val inputFolder = new File("C:/Users/Administrator/Desktop/MyScalismoProject/data/ArotaModels")
  val outputFolder = new File("C:/Users/Administrator/Desktop/MyScalismoProject/data/AlignedModels")
  outputFolder.mkdirs()

  val numSamplePoints = 50000  // 50000 points for high-quality aortic details
  val rng = new Random(42)

  def triangleArea(p1: Point[_3D], p2: Point[_3D], p3: Point[_3D]): Double = {
    val a = p2 - p1
    val b = p3 - p1
    val cross = EuclideanVector3D(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x
    )
    0.5 * cross.norm
  }

  // Improved sampling method: area-based stratified sampling
  def sample(mesh: TriangleMesh3D, n: Int): IndexedSeq[Point[_3D]] = {
    println(s"  Starting to sample ${n} points...")
    
    // Calculate area and cumulative area for each triangle
    val triangles = mesh.triangulation.triangleIds.map(tid => {
      val tri = mesh.triangulation.triangle(tid)
      val p1 = mesh.pointSet.point(tri.ptId1)
      val p2 = mesh.pointSet.point(tri.ptId2)
      val p3 = mesh.pointSet.point(tri.ptId3)
      (tid, triangleArea(p1, p2, p3), p1, p2, p3)
    }).toArray

    val totalArea = triangles.map(_._2).sum
    println(f"  Total mesh area: ${totalArea}%.2f")
    
    // Assign sampling point count for each triangle (based on area proportion)
    val trianglePointCounts = triangles.map { case (tid, area, _, _, _) =>
      val ratio = area / totalArea
      val pointCount = math.round(n * ratio).toInt  // Assign points based on area proportion
      (tid, pointCount, area)
    }
    
    val totalAssignedPoints = trianglePointCounts.map(_._2).sum
    println(s"  Initial total assigned points: ${totalAssignedPoints}")
    
    // Adjust point allocation to ensure total equals target point count
    val adjustedCounts = if (totalAssignedPoints != n) {
      val diff = n - totalAssignedPoints
      // Sort by area size, prioritize adjusting triangles with larger areas
      val sortedIndices = triangles.zipWithIndex.sortBy(-_._1._2).map(_._2)
      
      trianglePointCounts.zipWithIndex.map { case ((tid, count, area), idx) =>
        if (diff > 0 && sortedIndices.take(diff).contains(idx)) {
          (tid, count + 1, area)  // Add 1 point
        } else if (diff < 0 && count > 0 && sortedIndices.take(-diff).contains(idx)) {
          (tid, math.max(0, count - 1), area)  // Reduce 1 point, but not less than 0
        } else {
          (tid, count, area)
        }
      }
    } else trianglePointCounts
    
    val finalTotalPoints = adjustedCounts.map(_._2).sum
    println(s"  Final total assigned points: ${finalTotalPoints}")
    
    // Uniform sampling within each triangle
    def sampleInTriangleUniform(p1: Point[_3D], p2: Point[_3D], p3: Point[_3D], numPoints: Int): Seq[Point[_3D]] = {
      if (numPoints <= 0) return Seq.empty
      
      (0 until numPoints).map { i =>
        // Use more uniform sampling method
        val u = rng.nextDouble()
        val v = rng.nextDouble()
        
        // Ensure points are within triangle
        val (a, b) = if (u + v > 1.0) (1.0 - u, 1.0 - v) else (u, v)
        val c = 1.0 - a - b
        
        Point(
          a * p1.x + b * p2.x + c * p3.x,
          a * p1.y + b * p2.y + c * p3.y,
          a * p1.z + b * p2.z + c * p3.z
        )
      }
    }
    
    // Sample corresponding number of points from each triangle
    val sampledPoints = adjustedCounts.zipWithIndex.flatMap { case ((tid, pointCount, _), idx) =>
      val (_, _, p1, p2, p3) = triangles(idx)
      sampleInTriangleUniform(p1, p2, p3, pointCount)
    }
    
    println(s"  Actual sampled points: ${sampledPoints.length}")
    
    // Unify point order: sort by spatial coordinates (Z first, then Y, finally X)
    val sortedPoints = sampledPoints.sortBy(p => (p.z, p.y, p.x))
    println(s"  Points sorted by spatial coordinates")
    
    sortedPoints.toIndexedSeq
  }

  // Anatomical landmark-based alignment function
  def alignToAnatomicalLandmark(points: IndexedSeq[Point[_3D]]): IndexedSeq[Point[_3D]] = {
    // Find the estimated position of aortic root (heart connection point)
    // Usually the point with highest Z coordinate and in the central area of the model
    val zMax = points.map(_.z).max
    val zMin = points.map(_.z).min
    val zRange = zMax - zMin
    
    // Find points within the upper 25% Z coordinate range
    val upperZThreshold = zMax - zRange * 0.25
    val upperPoints = points.filter(_.z >= upperZThreshold)
    
    // Among these upper points, find the one closest to XY plane center as aortic root
    val xCenter = points.map(_.x).sum / points.length
    val yCenter = points.map(_.y).sum / points.length
    
    val aorticRoot = upperPoints.minBy { p =>
      val dx = p.x - xCenter
      val dy = p.y - yCenter
      math.sqrt(dx * dx + dy * dy)
    }
    
    println(f"  Detected aortic root position: (${aorticRoot.x}%.2f, ${aorticRoot.y}%.2f, ${aorticRoot.z}%.2f)")
    println(f"  Z range: ${zMin}%.2f to ${zMax}%.2f")
    
    // Move aortic root to origin
    points.map(p => Point(
      p.x - aorticRoot.x, 
      p.y - aorticRoot.y, 
      p.z - aorticRoot.z
    ))
  }

  // Add scale normalization function
  def normalizeScale(points: IndexedSeq[Point[_3D]]): IndexedSeq[Point[_3D]] = {
    val distances = points.map(p => math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z))
    val maxDistance = distances.max
    val avgDistance = distances.sum / distances.length
    
    println(f"  Max distance: ${maxDistance}%.2f, Average distance: ${avgDistance}%.2f")
    
    // Normalize to average distance of 1
    val scaleFactor = 1.0 / avgDistance
    points.map(p => Point(p.x * scaleFactor, p.y * scaleFactor, p.z * scaleFactor))
  }

  val files = inputFolder.listFiles().filter(_.getName.endsWith(".stl")).sorted
  
  println(s"Starting to process ${files.length} STL files with anatomical position alignment...")
  
  for (file <- files) {
    println(s"\nProcessing file: ${file.getName}")
    val mesh = MeshIO.readMesh(file).get
    val rawPoints = sample(mesh, numSamplePoints)
    
    // Apply anatomical alignment processing
    val alignedPoints = alignToAnatomicalLandmark(rawPoints) 
    val normalizedPoints = normalizeScale(alignedPoints)
    
    val csv = new File(outputFolder, s"anatomical_${file.getName.stripSuffix(".stl")}.csv")
    val writer = new PrintWriter(csv)
    normalizedPoints.foreach(p => writer.println(s"${p.x},${p.y},${p.z}"))
    writer.close()
    println(s"${csv.getName} saved (anatomically aligned and scale normalized)")
  }
  
  println(s"\nCompleted! All models have been anatomically aligned and scale normalized based on landmark points.")
}
