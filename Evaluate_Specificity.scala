//> using scala "3.3"
//> using dep "org.scalanlp::breeze:2.1.0"
//> using dep "org.scalanlp::breeze-natives:2.1.0"

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Paths}
import scala.io.Source
import scala.util.Random
import breeze.linalg.*
import breeze.numerics.*
import breeze.stats.*

object Evaluate_specificity extends App {

  // ====== 配置（按需修改） ======
  val alignedDir   = "C:/Users/Administrator/Desktop/MyScalismoProject/data/AlignedModels"
  val Ks           = (1 to 3).toVector       // 评估的K范围（测试用小范围）
  val Kmax         = Ks.max
  val S            = 10                       // 每个K采样的形状个数（测试用小数量）
  val maxPoints    = 15000                    // 评估用点数（与PCA阶段一致以避免内存问题）
  val seed         = 20250809
  val useMahalanobisConstraint = true         // 是否使用马氏半径约束
  val tau          = 2.0                      // 马氏半径阈值 d_M <= tau（常用2）
  val outDir       = "specificity_out"        // 输出目录

  new File(outDir).mkdirs()
  val rng = new Random(seed)

  // ====== 读取 CSV 并构造 X (M x D)，D=3*points ======
  val files = new File(alignedDir).listFiles().filter(_.getName.toLowerCase.endsWith(".csv")).sorted
  require(files.nonEmpty, s"No CSV found in $alignedDir")
  val M = files.length

  def countLines(path: String): Int = {
    val src = Source.fromFile(path)("UTF-8")
    try src.getLines().length finally src.close()
  }
  val minLines = files.map(f => countLines(f.getAbsolutePath)).min
  val take = math.min(maxPoints, minLines)
  val D = take * 3
  println(s"[Load] M=$M, points=$take, D=$D")

  def readXYZ(path: String, takePoints: Int): DenseVector[Double] = {
    val src = Source.fromFile(path)("UTF-8")
    val it = src.getLines()
    val buf = new Array[Double](takePoints * 3)
    var idx = 0
    var taken = 0
    try {
      while (it.hasNext && taken < takePoints) {
        val parts = it.next().split(",")
        if (parts.length >= 3) {
          buf(idx)   = parts(0).toDouble
          buf(idx+1) = parts(1).toDouble
          buf(idx+2) = parts(2).toDouble
          idx   += 3
          taken += 1
        }
      }
    } finally src.close()
    DenseVector(buf)
  }

  val X = DenseMatrix.zeros[Double](M, D)
  for (i <- 0 until M) {
    X(i, ::) := readXYZ(files(i).getAbsolutePath, take).t
    if (i % 10 == 0) println(s"[Load] $i / $M")
  }

  // 预先计算每个真实样本的平方范数，用于快速距离
  val rowNorm2 = DenseVector.zeros[Double](M)
  for (i <- 0 until M) {
    val row = X(i, ::).t
    rowNorm2(i) = norm(row, 2) * norm(row, 2)
  }

  // ====== PCA（Gram 技巧） ======
  // 中心化
  def colMean(mat: DenseMatrix[Double]): DenseVector[Double] = {
    val result = DenseVector.zeros[Double](mat.cols)
    for (j <- 0 until mat.cols) {
      var sum = 0.0
      for (i <- 0 until mat.rows) {
        sum += mat(i, j)
      }
      result(j) = sum / mat.rows
    }
    result
  }
  val mu = colMean(X)
  val Xc = X.copy
  for (row <- 0 until Xc.rows) {
    Xc(row, ::) := Xc(row, ::) - mu.t
  }

  // G = Xc * Xc^T   (M x M)，特征分解
  val G = Xc * Xc.t
  val eigResult = eigSym(G)
  val evalsAsc = eigResult.eigenvalues
  val evecs = eigResult.eigenvectors
  
  // 按特征值降序排列
  val evalArray = evalsAsc.toArray
  val indices = evalArray.zipWithIndex.sortBy(-_._1).map(_._2)
  val evals = DenseVector(indices.map(evalArray))
  val U = DenseMatrix.zeros[Double](evecs.rows, evecs.cols)
  for (i <- indices.indices) {
    U(::, i) := evecs(::, indices(i))
  }

  // 奇异值 σ_k = sqrt(λ_Gk)，其中 Cov = Xc^T Xc / (M-1)，其特征值 = σ_k^2 / (M-1)
  val sig = DenseVector(evals.toArray.map(e => math.sqrt(math.max(e, 1e-10))))
  
  // V_full = Xc^T * U * diag(1/sig)   => (D x M)，列正交单位
  val XcT_U = Xc.t * U
  val Vfull = DenseMatrix.zeros[Double](D, M)
  for (i <- 0 until M) {
    val scale = if (sig(i) > 1e-12) 1.0 / sig(i) else 0.0
    for (j <- 0 until D) {
      Vfull(j, i) = XcT_U(j, i) * scale
    }
  }

  def specificityForK(K: Int): (DenseVector[Double], Double, Double, Int, Int) = {
    val V = Vfull(::, 0 until K).toDenseMatrix        // D x K
    // 真实协方差的前K特征值：lambda_k = σ_k^2 / (M-1)
    val sigK = sig(0 until K)
    val lambdaK = DenseVector.zeros[Double](K)
    for (i <- 0 until K) {
      lambdaK(i) = sigK(i) * sigK(i) / (M - 1).toDouble
    }
    val sqrtLambdaK = DenseVector(lambdaK.toArray.map(math.sqrt))
    val out = DenseVector.zeros[Double](S)

    var accepted = 0
    var clipped  = 0

    for (s <- 0 until S) {
      var z = DenseVector.zeros[Double](K)
      var ok = false
      var tries = 0
      // 采样 z ~ N(0, I)，可选马氏约束 d_M = ||z|| <= tau
      while (!ok && tries < 50) {
        for (i <- 0 until K) {
          z(i) = rng.nextGaussian()
        }
        if (!useMahalanobisConstraint || norm(z, 2) <= tau) ok = true
        tries += 1
      }
      if (!ok && useMahalanobisConstraint) {
        // clip 到球面上
        val znorm = norm(z, 2)
        for (i <- 0 until K) {
          z(i) = z(i) / znorm * tau
        }
        clipped += 1
      } else {
        accepted += 1
      }

      val b = DenseVector.zeros[Double](K)
      for (i <- 0 until K) {
        b(i) = sqrtLambdaK(i) * z(i)
      }
      val x = mu + (V * b)                            // 生成形状（D维）
      val xNorm2 = norm(x, 2) * norm(x, 2)

      // 快速距离：||x - y||^2 = ||x||^2 + ||y||^2 - 2 x·y
      var minD2 = Double.MaxValue
      for (i <- 0 until M) {
        val dot = X(i, ::).t dot x
        val d2 = rowNorm2(i) + xNorm2 - 2.0 * dot
        if (d2 < minD2) minD2 = d2
      }
      val rmse = math.sqrt(math.max(0.0, minD2)) / math.sqrt(D.toDouble)
      out(s) = rmse
    }

    val mean = out.toArray.sum / out.length
    val variance = out.toArray.map(x => (x - mean) * (x - mean)).sum / (out.length - 1)
    val sd = math.sqrt(variance)
    (out, mean, sd, accepted, clipped)
  }

  // ====== 逐K评估并写出 ======
  val summary = new StringBuilder("K,S,mean,std,accepted,clipped,tau,mahalanobis\n")
  Ks.foreach { K =>
    println(s"[Specificity] K=$K ...")
    val (vals, mean, sd, acc, clip) = specificityForK(K)

    // 明细 CSV：每个采样的最近RMSE
    val detail = new StringBuilder("sample_id,nearest_rmse\n")
    for (i <- 0 until S) detail.append(s"$i,${vals(i)}\n")
    Files.write(Paths.get(s"$outDir/specificity_K${K}.csv"), detail.toString.getBytes("UTF-8"))

    summary.append(s"$K,$S,$mean,$sd,$acc,$clip,$tau,$useMahalanobisConstraint\n")
    println(f" -> mean=$mean%.6f  std=$sd%.6f  accepted=$acc  clipped=$clip")
  }
  Files.write(Paths.get(s"$outDir/specificity_summary.csv"), summary.toString.getBytes("UTF-8"))

  println(s"[Done] Results written to: $outDir")
}
