//> using scala "3.3"
//> using dep "org.scalanlp::breeze:2.1.0"
//> using dep "org.scalanlp::breeze-natives:2.1.0"

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Paths}
import scala.io.Source
import breeze.linalg.*
import breeze.numerics.*
import breeze.stats.*

/**
 * 统计形状模型泛化能力评估工具
 * 使用留一法(LOO)交叉验证评估不同主成分个数下的重建误差
 * 生成图8B所需的数据
 */
object Evaluate extends App {

  // ====== 配置参数 ======
  val alignedDir = "C:/Users/Administrator/Desktop/MyScalismoProject/src/main/scala/data/AlignedModels"  // 对齐后的点云CSV目录
  val outputFile = "evaluation_generalization_loormse.csv"  // 输出文件名
  val Ks         = (5 to 50 by 5).toVector  // 评估的主成分个数范围: 5, 10, 15, ..., 50
  val Kmax       = Ks.max  // 最大主成分个数
  val maxPoints  = 5000    // 最大点数限制，用于节省内存

  // ====== 读取 CSV 并构造 (M x D) 矩阵，D = 3 * maxPoints ======
  println("统计形状模型泛化能力评估")
  println(s"数据目录: $alignedDir")
  println(s"最大点数: $maxPoints")
  println(s"评估的K值范围: ${Ks.mkString(", ")}")
  println("============================================================")

  val files = new File(alignedDir).listFiles().filter(f => f.getName.toLowerCase.endsWith(".csv")).sorted
  require(files.nonEmpty, s"错误：未找到CSV文件！目录: $alignedDir")
  val M = files.length
  println(s"找到 $M 个CSV文件")

  def readXYZ(path: String, takePoints: Int): DenseVector[Double] = {
    // 假设每行是 "x,y,z"
    val it = Source.fromFile(path)("UTF-8").getLines()
    val buf = new Array[Double](takePoints * 3)
    var idx = 0
    var taken = 0
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
    DenseVector(buf) // 长度 = 3 * takePoints
  }

  // 先读第一份，确认可取的点数
  val firstCount = Source.fromFile(files.head)("UTF-8").getLines().length
  val take = math.min(maxPoints, firstCount)
  val D = take * 3

  val X = DenseMatrix.zeros[Double](M, D)
  for (i <- 0.until(M)) {
    val v = readXYZ(files(i).getAbsolutePath, take)
    X(i, ::) := v.t
    if (i % 10 == 0) println(s"[load] $i / $M")
  }
  println(s"Loaded matrix: M=$M, D=$D (points=$take)")

  // ====== 小工具 ======
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

  // 记录每个样本、每个K的RMSE
  val rmse = DenseMatrix.zeros[Double](M, Ks.length)

  // ====== 精确 LOO（Gram 技巧） ======
  for (i <- 0.until(M)) {
    // 1) 训练集（去掉第 i 行）
    val trainIdx = (0 until M).filter(_ != i)
    val Xtrain = DenseMatrix.zeros[Double](trainIdx.length, D)
    for ((row, idx) <- trainIdx.zipWithIndex) {
      Xtrain(idx, ::) := X(row, ::)
    }

    // 2) 均值 + 中心化
    val mu = colMean(Xtrain)
    val Xc = Xtrain.copy
    for (row <- 0 until Xc.rows) {
      Xc(row, ::) := Xc(row, ::) - mu.t
    }

    // 3) Gram：G = Xc * Xc^T  ( (M-1) x (M-1) )
    val G = Xc * Xc.t
    val eigResult = eigSym(G)
    val evals = eigResult.eigenvalues
    val evecs = eigResult.eigenvectors
    
    // 按特征值降序排列
    val evalArray = evals.toArray
    val indices = evalArray.zipWithIndex.sortBy(-_._1).map(_._2)
    val validK = math.min(Kmax, evalArray.count(_ > 1e-10))
    val top = indices.take(validK)
    
    // 手动构建选中的特征向量矩阵
    val U = DenseMatrix.zeros[Double](evecs.rows, validK)
    val selectedEvals = DenseVector.zeros[Double](validK)
    for (i <- top.indices) {
      val idx = top(i)
      selectedEvals(i) = math.max(evalArray(idx), 1e-10)
      U(::, i) := evecs(::, idx)
    }
    val sig = sqrt(selectedEvals)

    // 数值安全：避免除0
    val invSig = DenseVector(sig.toArray.map(s => if (s > 1e-12) 1.0/s else 0.0))

    // 4) PCs: V = Xc^T * U * diag(1/sig)   => (D x validK)
    val XcT_U = Xc.t * U
    val Vfull = DenseMatrix.zeros[Double](D, validK)
    for (i <- 0 until validK) {
      val scale = if (sig(i) > 1e-12) 1.0 / sig(i) else 0.0
      for (j <- 0 until D) {
        Vfull(j, i) = XcT_U(j, i) * scale
      }
    }

    // 5) 测试样本 i 重建/误差（逐K）
    val xTest = X(i, ::).t
    val xC    = xTest - mu
    for ((k, j) <- Ks.zipWithIndex) {
      val actualK = math.min(k, validK)
      if (actualK > 0) {
        // 手动构建V矩阵的子集
        val V = Vfull(::, 0 until actualK).toDenseMatrix
        val b = V.t * xC
        val reconstruction = V * b
        val xHat = mu + reconstruction
        val err = norm(xHat - xTest) / math.sqrt(D)
        rmse(i, j) = err
      } else {
        rmse(i, j) = Double.MaxValue
      }
    }
    println(s"[LOO] $i / $M done.")
  }

  // ====== 汇总 & 保存 ======
  val meanByK = mean(rmse(::, *)).t
  val stdByK  = stddev(rmse(::, *)).t
  println("========== 泛化能力评估结果 ==========")
  println("K, LOO_RMSE_mean, LOO_RMSE_std")
  Ks.indices.foreach { j => 
    println(f"${Ks(j)}, ${meanByK(j)}%.6f, ${stdByK(j)}%.6f")
  }

  // 将结果写入CSV文件
  val out = new StringBuilder("K,mean,std\n")
  Ks.indices.foreach { j => out.append(f"${Ks(j)},${meanByK(j)}%.8f,${stdByK(j)}%.8f\n") }
  Files.write(Paths.get(outputFile), out.toString.getBytes("UTF-8"))
  
  println(s"\n结果已保存到: $outputFile")
  println("评估完成！")
}
