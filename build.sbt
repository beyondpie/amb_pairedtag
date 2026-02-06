// check sbt build settings
// https://www.scala-sbt.org/1.x/docs/Basic-Def.html
// https://www.scala-sbt.org/1.x/docs/Task-Graph.html
// https://www.scala-sbt.org/1.x/docs/Scopes.html

// For parallel settings:
// https://www.scala-sbt.org/1.x/docs/Parallel-Execution.html
Global / semanticdbEnabled := true
Global / onChangedBuildSource := ReloadOnSourceChanges
// Build-wide settings
ThisBuild / organization := "io.github.beyondpie"
ThisBuild / organizationName := "zulab"
ThisBuild / scalaVersion := "3.7.3"
ThisBuild / logLevel := Level.Info
ThisBuild / resolvers += "Bioviz".at(
  "https://nexus.bioviz.org/repository/maven-releases/")

ThisBuild / Compile / scalacOptions := List(
  "-encoding",
  "utf8",
  "-feature",
  "-language:implicitConversions",
  "-language:existentials",
  "-unchecked",
  "-explain-types",
  "-explain",
  "-deprecation",
  "-preview"
)

// Ref: https://www.scala-sbt.org/1.x/docs/Multi-Project.html
lazy val bioscala = (project in file("bioscala"))
  .settings(
    name := "bioscala",
    version := "0.7.0",
    Test / logBuffered := false,
    javaOptions ++= Seq(
      "Dscala.concurrent.context.numThreads=1",
      "Dscala.concurrent.context.maxThreads=4"
    ),
    libraryDependencies ++= Seq(
      "org.scalatest" %% "scalatest" % "3.2.19" % "test",
      "com.lihaoyi" %% "os-lib" % "0.11.4",
      "org.scala-lang.modules" %% "scala-parallel-collections" % "1.2.0",
      "com.github.samtools" % "htsjdk" % "4.2.0",
      // only works on scala 3.6.3 now
      // "com.lihaoyi" % "ammonite" % "3.0.2" % "test" cross CrossVersion.full,
      "commons-io"           % "commons-io"  % "2.19.0",
      "org.jsoup"            % "jsoup"       % "1.20.1",
      // FIXME: our codes used uncompatible ones in 4.3.0
      "com.github.haifengl" %% "smile-scala" % "4.2.0",
      "org.bytedeco" % "javacpp" % "1.5.11",
      "org.bytedeco" % "openblas" % "0.3.28-1.5.11",
      "org.bytedeco" % "arpack-ng" % "3.9.1-1.5.11",
      "org.slf4j" % "slf4j-simple" % "2.0.17",
      "org.slf4j" % "slf4j-api"    % "2.0.17",
      "org.broad.igv" % "bigwig" % "3.0.0",
      "dev.optics" %% "monocle-core"  % "3.3.0",
      "dev.optics" %% "monocle-macro" % "3.3.0",
      "org.apache.commons" % "commons-math3" % "3.6.1"
    )
  )

lazy val pairedtag = (project in file("100.project"))
  .dependsOn(bioscala)
  .settings(
    name := "AdultMouseBrainPairedTag",
    version := "0.9"
  )
lazy val DataPreprocess = (project in file("00.datapreprocess"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "DataPrreprocess",
    version := "0.9"
  )

lazy val Integration = (project in file("03.integration"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "Integration",
    version := "1.0"
  )

lazy val CallPeak = (project in file("04.peakcalling"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "CallPeak",
    version := "1.0"
  )

lazy val CRE = (project in file("05.CRE"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "CRE",
    version := "1.0"
  )

lazy val ChromHMM = (project in file("06.ChromHMM"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "ChromHMM",
    version := "1.0"
  )


lazy val DL = (project in file("07.deeplearning"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "DL",
    version := "0.5"
  )


lazy val Epimem = (project in file("09.epimem"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "Epimem",
    version := "0.5"
  )


lazy val SuperEnhancer = (project in file("10.superEnhancer"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "SuperEnhancer",
    version := "1.0"
  )

lazy val DE = (project in file("12.DE"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "DE",
    version := "1.0",
    Test/fork := true,
    Test/javaOptions ++= Seq(
      "-Dscala.concurrent.context.numThreads1",
      "-Dscala.concurrent.context.maxThreads=2"
    )
  )


lazy val Cicero = (project in file("13.cicero"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "Cicero",
    version := "1.0"
  )

lazy val repressiveMarks = (project in file("17.repressiveMarks"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "repressiveMarks",
    version := "0.5.99"
  )


lazy val sexDE = (project in file("18.sexDE"))
  .dependsOn(bioscala, pairedtag)
  .settings(
    name := "sexDE",
    version := "0.1"
  )




