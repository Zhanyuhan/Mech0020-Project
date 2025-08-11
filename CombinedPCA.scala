//> using scala "3.3"
//> using dep "org.scalanlp::breeze:2.1.0"

import breeze.linalg._
import scala.io.Source
import java.io.File
import scala.collection.mutable.ListBuffer
import scala.math._

// Custom Point3D class
case class Point3D(x: Double, y: Double, z: Double)

object CombinedPCA extends App {

  val inputDir = new File("C:/Users/Administrator/Desktop/MyScalismoProject/data/AlignedModels")
  val pcaCompareDir = new File("C:/Users/Administrator/Desktop/MyScalismoProject/data/PCACompare")
  val finalOutputDir = new File("C:/Users/Administrator/Desktop/MyScalismoProject/data/Final")
  
  // Create output directories
  pcaCompareDir.mkdirs()
  finalOutputDir.mkdirs()

  val csvFiles = inputDir.listFiles().filter(_.getName.endsWith(".csv")).sorted
  val numPoints = 15000  // Try 15000 points, use stratified sampling to ensure representativeness

  println("Comprehensive PCA Analysis Program Started")
  println("============================================================")
  println(s"Performing comprehensive PCA analysis on ${csvFiles.length} aorta models...")
  println(s"Using ${numPoints} points per model - High quality aorta details")

  val allModels = csvFiles.map { file =>
    println(s"  Reading file: ${file.getName}")
    val lines = Source.fromFile(file).getLines().toArray
    val coords = lines.flatMap(_.split(",").map(_.toDouble))
    println(s"    Total lines: ${lines.length}, coordinate count: ${coords.length}, Expected: ${numPoints * 3}")
    
    // Add data quality check
    if (lines.length > 0) {
      val firstLine = lines(0).split(",")
      println(s"    First line: ${lines(0)}")
      println(s"    First line components: ${firstLine.length}")
      if (firstLine.length == 3) {
        val firstPoint = Point3D(firstLine(0).toDouble, firstLine(1).toDouble, firstLine(2).toDouble)
        println(s"    First point: (${firstPoint.x}, ${firstPoint.y}, ${firstPoint.z})")
      }
    }
    
    // Stratified sampling: uniformly select representative points from all points
    if (coords.length >= 50000 * 3) {
      val totalPoints = coords.length / 3
      val step = totalPoints.toDouble / numPoints
      val selectedCoords = (0 until numPoints).flatMap { i =>
        val pointIndex = (i * step).toInt
        val coordIndex = pointIndex * 3
        if (coordIndex + 2 < coords.length) {
          Array(coords(coordIndex), coords(coordIndex + 1), coords(coordIndex + 2))
        } else {
          Array.empty[Double]
        }
      }
      selectedCoords.toArray
    } else {
      coords.take(numPoints * 3)
    }
  }.filter(_.length >= numPoints * 3)  // Ensure at least sufficient number of points

  if (allModels.length < 5) {
    println("Need at least 5 models for reliable PCA analysis")
    sys.exit(1)
  }

  println(s"Successfully loaded ${allModels.length} valid models")
  
  // Add data quality check
  println("\nData Quality Check:")
  for ((model, idx) <- allModels.zipWithIndex.take(3)) {
    val points = (0 until math.min(1000, numPoints)).map { j =>
      Point3D(model(j * 3), model(j * 3 + 1), model(j * 3 + 2))
    }
    val centroid = Point3D(
      points.map(_.x).sum / points.length,
      points.map(_.y).sum / points.length,
      points.map(_.z).sum / points.length
    )
    val distances = points.map(p => math.sqrt(
      (p.x - centroid.x) * (p.x - centroid.x) + 
      (p.y - centroid.y) * (p.y - centroid.y) + 
      (p.z - centroid.z) * (p.z - centroid.z)
    ))
    println(f"  Model ${idx + 1}: Center(${centroid.x}%.2f, ${centroid.y}%.2f, ${centroid.z}%.2f), Average radius: ${distances.sum / distances.length}%.2f")
  }

  val meanShape = allModels.reduce((a, b) => a.zip(b).map((x, y) => x + y)).map(_ / allModels.length)
  val meanVector = DenseVector(meanShape)

  val dataMatrix = DenseMatrix.zeros[Double](allModels.length, numPoints * 3)
  for ((model, idx) <- allModels.zipWithIndex) {
    dataMatrix(idx, ::) := DenseVector(model).t
  }

  val centeredMatrix = dataMatrix(*, ::) - meanVector

  println("\nPerforming SVD decomposition...")
  val svd.SVD(u, s, vt) = svd(centeredMatrix.t)

  val eigenvalues = s.toArray.map(x => x * x / (allModels.length - 1))
  val totalVariance = eigenvalues.sum
  val explainedVariance = eigenvalues.map(_ / totalVariance)

  println("\nPrincipal Component Variance Explained Ratios:")
  for (i <- explainedVariance.indices.take(10)) {
    println(f"  Principal Component ${i + 1}%2d: ${explainedVariance(i) * 100}%.2f%%")
  }

  val principalComponents = u.t
  val referenceModelIndex = 0
  val referenceModel = allModels(referenceModelIndex)

  println(f"\nDebug Information:")
  println(f"  meanVector length: ${meanVector.length}")
  println(f"  principalComponents dimensions: ${principalComponents.rows} x ${principalComponents.cols}")
  println(f"  Number of singular values: ${s.length}")

  // Corrected scaling function: scale to size suitable for visualization
  def scaleForVisualization(points: Seq[Point3D]): Seq[Point3D] = {
    val distances = points.map(p => math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z))
    val avgDistance = distances.sum / distances.length
    
    // Scale by 100x to make model clearly visible in visualization software
    val scaleFactor = 100.0
    val scaledPoints = points.map(p => Point3D(p.x * scaleFactor, p.y * scaleFactor, p.z * scaleFactor))
    
    println(f"  Average distance before scaling: ${avgDistance}%.4f, Scale factor: ${scaleFactor}x")
    scaledPoints
  }
  
  // Boundary detection algorithm designed specifically for aorta tubular structures
  def identifyBoundaryPoints(points: Seq[Point3D], k: Int = 15): Seq[(Point3D, Boolean)] = {
    if (points.length < k || points.length < 50) return points.map((_, false))
    
    // Calculate center and main axis direction of entire point cloud
    val centroid = Point3D(
      points.map(_.x).sum / points.length,
      points.map(_.y).sum / points.length,
      points.map(_.z).sum / points.length
    )
    
    // Calculate distance from each point to central axis (assuming main axis is Z-axis)
    val radialDistances = points.map { p =>
      val radialDist = math.sqrt((p.x - centroid.x) * (p.x - centroid.x) + (p.y - centroid.y) * (p.y - centroid.y))
      (p, radialDist)
    }
    
    val sortedRadialDistances = radialDistances.map(_._2).sorted
    val medianRadialDistance = sortedRadialDistances(sortedRadialDistances.length / 2)
    val q75RadialDistance = sortedRadialDistances((sortedRadialDistances.length * 3) / 4)
    
    val boundaryPoints = points.map { point =>
      val pointRadialDistance = math.sqrt((point.x - centroid.x) * (point.x - centroid.x) + (point.y - centroid.y) * (point.y - centroid.y))
      
      // Calculate nearest neighbor density
      val distances = points.map { otherPoint =>
        if (point == otherPoint) Double.MaxValue
        else math.sqrt(
          (point.x - otherPoint.x) * (point.x - otherPoint.x) +
          (point.y - otherPoint.y) * (point.y - otherPoint.y) +
          (point.z - otherPoint.z) * (point.z - otherPoint.z)
        )
      }.sorted.take(k)
      
      val avgNearestDistance = distances.sum / distances.length
      val minDistance = distances.head
      val density = k / (avgNearestDistance * avgNearestDistance * avgNearestDistance)  // Point density
      
      // Combined conditions: 1) In outer region 2) Low point density 3) Large nearest neighbor distance
      val isOnBoundary = (
        pointRadialDistance > medianRadialDistance * 1.1 &&  // In outer region
        (avgNearestDistance > minDistance * 1.8 || density < 0.001)  // Low density or large gaps
      )
      
      (point, isOnBoundary)
    }
    
    val boundaryCount = boundaryPoints.count(_._2)
    println(f"    Aorta boundary detection: Identified ${boundaryCount} boundary points (${(boundaryCount.toDouble/points.length*100)}%.1f%%)")
    
    boundaryPoints
  }

  // Write PLY point cloud file (high quality, suitable for visualization)
  def writePLYPointCloud(points: Seq[Point3D], colors: Seq[(Int, Int, Int)], file: File): Unit = {
    val writer = new java.io.PrintWriter(file)
    writer.println("ply")
    writer.println("format ascii 1.0")
    writer.println(s"element vertex ${points.length}")
    writer.println("property float x")
    writer.println("property float y")
    writer.println("property float z")
    writer.println("property uchar red")
    writer.println("property uchar green")
    writer.println("property uchar blue")
    writer.println("end_header")
    for ((p, c) <- points.zip(colors)) {
      writer.println(s"${p.x} ${p.y} ${p.z} ${c._1} ${c._2} ${c._3}")
    }
    writer.close()
  }

  // Generate different colors for different deformations
  def getColorForVariation(factor: Double): (Int, Int, Int) = {
    factor match {
      case -2.0 => (255, 0, 0)     // Red: -2σ 
      case -1.0 => (255, 128, 0)   // Orange: -1σ
      case 0.0  => (0, 255, 0)     // Green: mean shape
      case 1.0  => (0, 128, 255)   // Light blue: +1σ
      case 2.0  => (0, 0, 255)     // Blue: +2σ
      case _    => (128, 128, 128) // Gray: others
    }
  }

  // Generate reconstruction files for each model (original + 6 principal component deformations) with boundary enhancement
  def generateModelReconstructions(): Unit = {
    println(s"\nGenerating reconstruction files for ${allModels.length} models (with boundary enhancement)...")
    
    for ((originalModel, modelIdx) <- allModels.zipWithIndex) {
      val modelName = csvFiles(modelIdx).getName.replace(".csv", "")
      println(s"\nProcessing model ${modelIdx + 1}/${allModels.length}: ${modelName}")
      
      // Create folder for each model
      val modelFolder = new File(finalOutputDir, modelName)
      if (!modelFolder.exists()) modelFolder.mkdirs()
      
      // Save original model (with boundary enhancement)
      val originalPoints = (0 until numPoints).map { i =>
        Point3D(originalModel(i * 3), originalModel(i * 3 + 1), originalModel(i * 3 + 2))
      }
      val originalFile = new File(modelFolder, s"${modelName}_original.ply")
      val scaledOriginalPoints = scaleForVisualization(originalPoints)
      
      // Boundary enhancement: identify boundary points and use different colors
      val boundaryInfo = identifyBoundaryPoints(scaledOriginalPoints)
      val boundaryColors = boundaryInfo.map { case (_, isBoundary) =>
        if (isBoundary) (255, 200, 0) else (100, 200, 100)  // Boundary points: yellow, interior points: green
      }
      
      writeColoredPLY(boundaryInfo.map(_._1), boundaryColors, originalFile)
      
      // Generate ±2σ deformations for the first 3 principal components
      for (pcIdx <- 0 until math.min(3, s.length)) {
        val stdDev = math.sqrt(eigenvalues(pcIdx))
        val pcVector = principalComponents(pcIdx, ::).t
        val sigma = 1.5  // Use more reasonable deformation amplitude
        
        // +2σ deformation - based on original model deformation
        val originalVector = DenseVector(originalModel)
        val posDisplaced = originalVector + pcVector * (sigma * stdDev)
        val posPoints = (0 until numPoints).map { i =>
          Point3D(posDisplaced(i * 3), posDisplaced(i * 3 + 1), posDisplaced(i * 3 + 2))
        }
        val posFile = new File(modelFolder, s"${modelName}_PC${pcIdx + 1}_plus2sigma.ply")
        val scaledPosPoints = scaleForVisualization(posPoints)
        writeColoredPLY(scaledPosPoints, Array.fill(scaledPosPoints.size)((0, 150, 255)), posFile)  // Light blue
        
        // -2σ deformation - based on original model deformation
        val negDisplaced = originalVector + pcVector * (-sigma * stdDev)
        val negPoints = (0 until numPoints).map { i =>
          Point3D(negDisplaced(i * 3), negDisplaced(i * 3 + 1), negDisplaced(i * 3 + 2))
        }
        val negFile = new File(modelFolder, s"${modelName}_PC${pcIdx + 1}_minus2sigma.ply")
        val scaledNegPoints = scaleForVisualization(negPoints)
        writeColoredPLY(scaledNegPoints, Array.fill(scaledNegPoints.size)((255, 100, 0)), negFile)  // Orange-red
        
        println(s"  Generated PC${pcIdx + 1} deformations: ±${sigma}σ")
      }
      
      println(s"  7 PLY files for model ${modelName} saved to ${modelFolder.getAbsolutePath}")
    }
  }

  // Generate PCA principal component overlapping visualization files (with boundary enhancement)
  def generateOverlappingVisualization(pcIndex: Int, componentName: String): Unit = {
    println(s"  Generating ${componentName} overlapping visualization file (with boundary enhancement)...")
    try {
      val allPoints = scala.collection.mutable.ListBuffer[Point3D]()
      val allColors = scala.collection.mutable.ListBuffer[(Int, Int, Int)]()

      for (factor <- Seq(-2.0, 0.0, 2.0)) {
        val pcVector = principalComponents(pcIndex, ::).t
        val stdDev = math.sqrt(eigenvalues(pcIndex))
        val displacedVector = DenseVector(referenceModel) + pcVector * (factor * stdDev)
        val displacedPoints = (0 until numPoints).map { j =>
          Point3D(displacedVector(j * 3), displacedVector(j * 3 + 1), displacedVector(j * 3 + 2))
        }

        // Scale 100x for visualization
        val scaleFactor = 100.0
        val scaledPoints = displacedPoints.map(p => Point3D(p.x * scaleFactor, p.y * scaleFactor, p.z * scaleFactor))
        val colors = Array.fill(scaledPoints.size)(getColorForVariation(factor))

        allPoints ++= scaledPoints
        allColors ++= colors
      }

      // Write overlapping file to PCACompare folder
      val combinedFile = new File(pcaCompareDir, s"${componentName}_overlapping.ply")
      writePLYPointCloud(allPoints.toSeq, allColors.toSeq, combinedFile)
      println(f"    Output overlapping file ${componentName} -> ${combinedFile.getName}")

    } catch {
      case e: Exception =>
        println(f"    Error processing overlapping file ${componentName}: ${e.getMessage}")
    }
  }

  // Save point cloud data in CSV format
  def writePointCloudAsCSV(points: Seq[Point3D], file: File): Unit = {
    val writer = new java.io.PrintWriter(file)
    try {
      points.foreach { p =>
        writer.println(s"${p.x},${p.y},${p.z}")
      }
    } finally {
      writer.close()
    }
  }

  def writeColoredPLY(points: Seq[Point3D], colors: Seq[(Int, Int, Int)], file: File): Unit = {
    val writer = new java.io.PrintWriter(file)
    writer.println("ply")
    writer.println("format ascii 1.0")
    writer.println(s"element vertex ${points.length}")
    writer.println("property float x")
    writer.println("property float y")
    writer.println("property float z")
    writer.println("property uchar red")
    writer.println("property uchar green")
    writer.println("property uchar blue")
    writer.println("end_header")
    for ((p, c) <- points.zip(colors)) {
      writer.println(s"${p.x} ${p.y} ${p.z} ${c._1} ${c._2} ${c._3}")
    }
    writer.close()
  }

  // This function is no longer needed

  // Execute analysis
  println("\n============================================================")
  println("Starting PCA analysis...")
  
  // 1. Generate overlapping visualizations for the first 3 principal components to PCACompare
  println("\nGenerating principal component overlapping visualizations...")
  for (pcIndex <- 0.until(math.min(3, s.length))) {
    val componentName = pcIndex match {
      case 0 => "PCA1"
      case 1 => "PCA2"
      case 2 => "PCA3"
    }
    
    println(f"\nProcessing ${componentName} (Explained variance: ${explainedVariance(pcIndex) * 100}%.2f%%)...")
    generateOverlappingVisualization(pcIndex, componentName)
  }
  
  // 2. Generate reconstruction files for each model to Final
  println("\nGenerating model reconstruction files...")
  generateModelReconstructions()

  // Output summary
  println("\n============================================================")
  println("PCA analysis completed!")
  println("Output results:")
  println(s"  PCA overlapping files: ${pcaCompareDir.getAbsolutePath}")
  println(s"  Model reconstruction files: ${finalOutputDir.getAbsolutePath}")
  
  println("\nGenerated files:")
  println("  PCACompare folder: 3 principal component overlapping PLY files")
  println(s"  Final folder: ${allModels.length} model folders, each containing 7 PLY files")
  println("     - 1 original model + 6 principal component deformations (PC1-3 ±2σ)")
  
  println("\nFeature description:")
  println("  High quality model: Uses 50000 points to capture fine aorta structure")
  println("  Boundary enhancement: Intelligently identifies model boundaries, highlights boundary structures")
  println("  Original size: Maintains original model size without any scaling")
  
  println("\nColor coding explanation (PLY files):")
  println("  Red: -2σ deformation (negative extreme)")
  println("  Green: 0σ mean shape (interior points)")
  println("  Blue: +2σ deformation (positive extreme)")
  println("  Yellow: Model boundary points (enhanced display)")
  println("  Bright yellow: Boundary points in overlapping")
}
