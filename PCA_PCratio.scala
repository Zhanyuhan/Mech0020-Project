//> using scala "3.3"
//> using dep "org.scalanlp:breeze_2.13:2.1.0"
//> using dep "org.scalanlp:breeze-natives_2.13:2.1.0"
//> using jvm "17"
//> using javaOpt "-Xmx16G"
//> using javaOpt "-Xms8G"
//> using javaOpt "-XX:+UseG1GC"
//> using javaOpt "-XX:G1HeapRegionSize=32M"
//> using javaOpt "-XX:MaxDirectMemorySize=4G"

import java.io.File
import breeze.linalg._
import breeze.numerics._
import scala.io.Source
import scala.util.Try

object PCAEigenvalues extends App {

  // 系统配置检查
  val runtime = Runtime.getRuntime
  val maxMemory = runtime.maxMemory() / (1024 * 1024 * 1024.0)
  val totalMemory = runtime.totalMemory() / (1024 * 1024 * 1024.0)
  val freeMemory = runtime.freeMemory() / (1024 * 1024 * 1024.0)
  
  println(f"🔧 系统内存配置:")
  println(f"   最大可用内存: ${maxMemory}%.2f GB")
  println(f"   当前分配内存: ${totalMemory}%.2f GB")
  println(f"   空闲内存: ${freeMemory}%.2f GB")
  
  // 路径和 CombinedPCA.scala 保持一致
  val dataDir = new File("C:/Users/Administrator/Desktop/MyScalismoProject/data/AlignedModels")
  val csvFiles = dataDir.listFiles().filter(_.getName.toLowerCase.endsWith(".csv")).sorted // 处理所有文件

  println(s"\n📂 准备处理所有 ${csvFiles.length} 个CSV文件...")
  
  // 智能采样策略 - 提高采样率以获得更准确的结果
  val availableMemoryGB = freeMemory + (maxMemory - totalMemory)
  val samplingStep = if (availableMemoryGB > 12) {
    5   // 高内存：每5个点采样1个 (约10000点，20%采样率)
  } else if (availableMemoryGB > 8) {
    10  // 中等内存：每10个点采样1个 (5000点，10%采样率)
  } else {
    20  // 低内存：每20个点采样1个 (2500点，5%采样率)
  }
  
  println(f"🧠 智能采样配置: 可用内存 ${availableMemoryGB}%.2f GB -> 采样步长 $samplingStep")
  
  // 分块读取策略，减少内存峰值
  val batchSize = 20  // 每批处理20个文件
  val batches = csvFiles.grouped(batchSize).toList
  
  println(s"📊 分批处理策略: ${batches.length}批，每批${batchSize}个文件")
  
  val allSampledData = batches.zipWithIndex.flatMap { case (batch, batchIdx) =>
    println(s"\n🔄 处理第${batchIdx + 1}/${batches.length}批文件...")
    
    batch.zipWithIndex.map { case (file, idx) =>
      val globalIdx = batchIdx * batchSize + idx + 1
      println(s"  读取文件 $globalIdx/${csvFiles.length}: ${file.getName}")
      
      val allValues = Source.fromFile(file).getLines().flatMap(_.split(",")).map(_.toDouble).toArray
      
      // 分层采样：保证采样点在整个模型上均匀分布
      val originalPoints = allValues.length / 3
      val sampledValues = (0 until originalPoints by samplingStep).flatMap { pointIdx =>
        val baseIdx = pointIdx * 3
        if (baseIdx + 2 < allValues.length) {
          Array(allValues(baseIdx), allValues(baseIdx + 1), allValues(baseIdx + 2))
        } else {
          Array.empty[Double]
        }
      }.toArray
      
      val sampledPoints = sampledValues.length / 3
      println(s"    原始点数: $originalPoints, 采样后: $sampledPoints (${(sampledPoints.toDouble/originalPoints*100).toInt}%)")
      sampledValues
    }
  }.toArray
  
  println(s"\n💾 构建数据矩阵...")
  val dataMatrix = DenseMatrix(allSampledData : _*)

  println(s"读取到 ${dataMatrix.rows} 个模型，每个模型 ${dataMatrix.cols} 个特征")
  println(s"数据矩阵大小: ${dataMatrix.rows} x ${dataMatrix.cols}")

  // 去均值
  println("计算均值...")
  val meanVector = breeze.stats.mean(dataMatrix(::, *))
  println("去中心化...")
  val centered = dataMatrix(*, ::) - meanVector.t

  println("执行SVD分解...")
  // 使用SVD方法，更稳定且不需要显式计算协方差矩阵
  val svd.SVD(u, s, vt) = breeze.linalg.svd(centered.t)
  
  // 计算特征值（奇异值的平方除以n-1）
  val eigenvalues = s.toArray.map(x => (x * x) / (dataMatrix.rows - 1).toDouble).sorted.reverse
  
  println(s"计算得到 ${eigenvalues.length} 个特征值")
  
  // 输出特征值及方差占比
  val totalVariance = eigenvalues.sum
  println("\n主成分特征值及方差占比：")
  println("=".repeat(70))
  
  var cumulativeVariance = 0.0
  eigenvalues.zipWithIndex.take(15).foreach { case (value, idx) =>
    val explained = value / totalVariance * 100
    cumulativeVariance += explained
    println(f"PC${idx + 1}%2d: 特征值 = $value%8.6f, 方差占比 = $explained%6.2f%%, 累计 = $cumulativeVariance%6.2f%%")
  }
  
  println("\n重要统计信息：")
  println("=".repeat(50))
  println(f"数据集规模: ${csvFiles.length} 个模型, 每个模型 ${dataMatrix.cols} 个特征")
  println(f"前3个主成分累计方差占比: ${eigenvalues.take(3).sum / totalVariance * 100}%.2f%%")
  println(f"前5个主成分累计方差占比: ${eigenvalues.take(5).sum / totalVariance * 100}%.2f%%")
  println(f"前10个主成分累计方差占比: ${eigenvalues.take(10).sum / totalVariance * 100}%.2f%%")
  
  // 保存结果到文件
  val outputFile = new File("PCA_eigenvalues_analysis.txt")
  val writer = new java.io.PrintWriter(outputFile)
  try {
    writer.println("PCA特征值分析报告")
    writer.println("=".repeat(50))
    writer.println(f"分析日期: ${java.time.LocalDateTime.now()}")
    writer.println(f"数据集规模: ${csvFiles.length} 个模型")
    writer.println(f"每个模型特征数: ${dataMatrix.cols}")
    writer.println(f"采样率: 每${samplingStep}个点采样1个")
    writer.println()
    
    writer.println("主成分特征值分析:")
    writer.println("-".repeat(70))
    
    var cumVar = 0.0
    eigenvalues.zipWithIndex.foreach { case (value, idx) =>
      val explained = value / totalVariance * 100
      cumVar += explained
      writer.println(f"PC${idx + 1}%2d: 特征值 = $value%10.6f, 方差占比 = $explained%6.2f%%, 累计 = $cumVar%6.2f%%")
    }
    
    writer.println()
    writer.println("重要统计:")
    writer.println(f"前3个主成分累计方差占比: ${eigenvalues.take(3).sum / totalVariance * 100}%.2f%%")
    writer.println(f"前5个主成分累计方差占比: ${eigenvalues.take(5).sum / totalVariance * 100}%.2f%%")
    writer.println(f"前10个主成分累计方差占比: ${eigenvalues.take(10).sum / totalVariance * 100}%.2f%%")
    
  } finally {
    writer.close()
  }
  
  println(f"\n分析结果已保存到: ${outputFile.getAbsolutePath}")
}
