ThisBuild / scalaVersion := "3.3.0"

// JVM内存设置
ThisBuild / javaOptions ++= Seq(
  "-Xmx20G",            // 最大堆内存12GB（你也可以用 8G，如果内存不够就提高）
  "-Xms2G",             // 初始堆内存2GB
  "-XX:+UseG1GC"        // 使用G1垃圾收集器
)

// 库依赖
libraryDependencies ++= Seq(
  "ch.unibas.cs.gravis" %% "scalismo-ui" % "0.92.0",
  "ch.unibas.cs.gravis" %% "scalismo"    % "0.92",
  "org.scalanlp"        %% "breeze"      % "2.1.0"
)
