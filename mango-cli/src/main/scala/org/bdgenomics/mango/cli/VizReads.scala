/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.mango.cli

import org.apache.hadoop.mapreduce.Job
import org.apache.spark.SparkContext
import org.bdgenomics.adam.models.VariantContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.cli._
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.adam.projections.AlignmentRecordField._
import org.bdgenomics.adam.projections.Projection
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.mango.models.OrderedTrackedLayout
import org.bdgenomics.formats.avro.{ AlignmentRecord, Feature, Genotype, GenotypeAllele, NucleotideContigFragment }
import org.fusesource.scalate.TemplateEngine
import org.json4s._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }
import org.scalatra.json._
import org.scalatra.ScalatraServlet

object VizReads extends ADAMCommandCompanion {
  val commandName: String = "viz"
  val commandDescription: String = "Genomic visualization for ADAM"
  var readsRefName = ""
  var variantsRefName = ""
  var reads: RDD[AlignmentRecord] = null
  var variants: RDD[Genotype] = null
  var reference: RDD[NucleotideContigFragment] = null
  var features: RDD[Feature] = null
  val trackHeight = 10
  val width = 1200
  val height = 400
  val base = 50

  def apply(cmdLine: Array[String]): ADAMCommand = {
    new VizReads(Args4j[VizReadsArgs](cmdLine))
  }

  def printTrackJson(layout: OrderedTrackedLayout[AlignmentRecord]): List[TrackJson] = {
    var tracks = new scala.collection.mutable.ListBuffer[TrackJson]
    for ((rec, track) <- layout.trackAssignments) {
      val aRec = rec._2.asInstanceOf[AlignmentRecord]
      tracks += new TrackJson(aRec.getReadName, aRec.getStart, aRec.getEnd, track)
    }
    tracks.toList
  }

  def printVariationJson(layout: OrderedTrackedLayout[Genotype]): List[VariationJson] = {
    var tracks = new scala.collection.mutable.ListBuffer[VariationJson]
    for ((rec, track) <- layout.trackAssignments) {
      val aRec = rec._2.asInstanceOf[Genotype]
      val referenceAllele = aRec.getVariant.getReferenceAllele
      val alternateAllele = aRec.getVariant.getAlternateAllele
      tracks += new VariationJson(aRec.getVariant.getContig.getContigName, aRec.getAlleles.mkString(" / "), aRec.getVariant.getStart, aRec.getVariant.getEnd, track)
    }
    tracks.toList
  }

  def printJsonFreq(array: Array[AlignmentRecord], region: ReferenceRegion): List[FreqJson] = {
    val freqMap = new java.util.TreeMap[Long, Long]
    var freqBuffer = new scala.collection.mutable.ListBuffer[FreqJson]

    // initiates map with 0 values
    var i0 = region.start
    while (i0 <= region.end) {
      freqMap.put(i0, 0)
      i0 += 1
    }

    // creates a point for each base showing its frequency
    for (rec <- array) {
      val aRec = rec.asInstanceOf[AlignmentRecord]
      var i = aRec.getStart
      while (i <= aRec.getEnd) {
        if (i >= region.start && i <= region.end)
          freqMap.put(i, freqMap(i) + 1)
        i += 1
      }
    }

    // convert to list of FreqJsons
    val iter = freqMap.keySet.iterator
    var key = 0L
    while (iter.hasNext) {
      key = iter.next()
      freqBuffer += FreqJson(key, freqMap(key))
    }

    freqBuffer.toList
  }
}

case class TrackJson(readName: String, start: Long, end: Long, track: Long)
case class VariationJson(contigName: String, alleles: String, start: Long, end: Long, track: Long)
case class FreqJson(base: Long, freq: Long)

class VizReadsArgs extends Args4jBase with ParquetArgs {
  @Argument(required = true, metaVar = "READS", usage = "The reads file to view", index = 0)
  var readPath: String = null

  @Argument(required = true, metaVar = "READS REFNAME", usage = "The reads reference to view", index = 1)
  var readsRefName: String = null

  @Argument(required = true, metaVar = "VARIANTS", usage = "The variants file to view", index = 2)
  var variantsPath: String = null

  @Argument(required = true, metaVar = "VARIANTS REFNAME", usage = "The variants reference to view", index = 3)
  var variantsRefName: String = null

  @Args4jOption(required = false, name = "-port", usage = "The port to bind to for visualization. The default is 8080.")
  var port: Int = 8080
}

class VizServlet extends ScalatraServlet with JacksonJsonSupport {
  protected implicit val jsonFormats: Formats = DefaultFormats
  var readsRegion = ReferenceRegion(VizReads.readsRefName, 0, 100)
  var variantsRegion = ReferenceRegion(VizReads.variantsRefName, 0, 100)

  get("/?") {
    redirect(url("reads"))
  }

  get("/reads") {
    contentType = "text/html"
    val readsRDD: RDD[AlignmentRecord] = VizReads.reads.filterByOverlappingRegion(readsRegion)
    val trackinput: RDD[(ReferenceRegion, AlignmentRecord)] = readsRDD.keyBy(ReferenceRegion(_))
    val filteredLayout = new OrderedTrackedLayout(trackinput.collect())
    val templateEngine = new TemplateEngine
    templateEngine.layout("mango-cli/src/main/webapp/WEB-INF/layouts/reads.ssp",
      Map("readsRegion" -> (readsRegion.referenceName, readsRegion.start.toString, readsRegion.end.toString),
        "width" -> VizReads.width.toString,
        "base" -> VizReads.base.toString,
        "numTracks" -> filteredLayout.numTracks.toString,
        "trackHeight" -> VizReads.trackHeight.toString))
  }

  get("/reads/:ref") {
    contentType = formats("json")
    val readRefRegion = ReferenceRegion(params("ref"), params("start").toLong, params("end").toLong)
    val readsRDD: RDD[AlignmentRecord] = VizReads.reads.filterByOverlappingRegion(readRefRegion)
    val trackinput: RDD[(ReferenceRegion, AlignmentRecord)] = readsRDD.keyBy(ReferenceRegion(_))
    val filteredLayout = new OrderedTrackedLayout(trackinput.collect())
    VizReads.printTrackJson(filteredLayout)
  }

  get("/overall") {
    contentType = "text/html"
    val readsRDD: RDD[AlignmentRecord] = VizReads.reads.filterByOverlappingRegion(readsRegion)
    val trackinput: RDD[(ReferenceRegion, AlignmentRecord)] = readsRDD.keyBy(ReferenceRegion(_))
    val filteredLayout = new OrderedTrackedLayout(trackinput.collect())
    val templateEngine = new TemplateEngine
    templateEngine.layout("mango-cli/src/main/webapp/WEB-INF/layouts/overall.ssp",
      Map("readsRegion" -> (readsRegion.referenceName, readsRegion.start.toString, readsRegion.end.toString),
        "variantsRegion" -> (variantsRegion.referenceName, variantsRegion.start.toString, variantsRegion.end.toString),
        "width" -> VizReads.width.toString,
        "base" -> VizReads.base.toString,
        "numTracks" -> filteredLayout.numTracks.toString,
        "trackHeight" -> VizReads.trackHeight.toString))
  }

  get("/freq") {
    contentType = "text/html"
    val templateEngine = new TemplateEngine
    templateEngine.layout("mango-cli/src/main/webapp/WEB-INF/layouts/freq.ssp",
      Map("readsRegion" -> (readsRegion.referenceName, readsRegion.start.toString, readsRegion.end.toString),
        "width" -> VizReads.width.toString,
        "height" -> VizReads.height.toString,
        "base" -> VizReads.base.toString))
  }

  get("/freq/:ref") {
    contentType = formats("json")
    val readRefRegion = ReferenceRegion(params("ref"), params("start").toLong, params("end").toLong)
    val filteredArray = VizReads.reads.filterByOverlappingRegion(readRefRegion).collect()
    VizReads.printJsonFreq(filteredArray, readRefRegion)
  }

  get("/variants") {
    contentType = "text/html"
    val variantsRDD: RDD[Genotype] = VizReads.variants.filterByOverlappingRegion(variantsRegion)
    val trackinput: RDD[(ReferenceRegion, Genotype)] = variantsRDD.keyBy(v => ReferenceRegion(ReferencePosition(v)))
    val filteredGenotypeTrack = new OrderedTrackedLayout(trackinput.collect())
    val templateEngine = new TemplateEngine
    val displayMap = Map("variantsRegion" -> (variantsRegion.referenceName, variantsRegion.start.toString, variantsRegion.end.toString),
      "width" -> VizReads.width.toString,
      "base" -> VizReads.base.toString,
      "numTracks" -> filteredGenotypeTrack.numTracks.toString,
      "trackHeight" -> VizReads.trackHeight.toString)
    templateEngine.layout("mango-cli/src/main/webapp/WEB-INF/layouts/variants.ssp",
      displayMap)
  }

  get("/variants/:ref") {
    contentType = formats("json")
    val variantRefRegion = ReferenceRegion(params("ref"), params("start").toLong, params("end").toLong)
    val variantsRDD: RDD[Genotype] = VizReads.variants.filterByOverlappingRegion(variantRefRegion)
    val trackinput: RDD[(ReferenceRegion, Genotype)] = variantsRDD.keyBy(v => ReferenceRegion(ReferencePosition(v)))
    val filteredGenotypeTrack = new OrderedTrackedLayout(trackinput.collect())
    VizReads.printVariationJson(filteredGenotypeTrack)
  }

}

class VizReads(protected val args: VizReadsArgs) extends ADAMSparkCommand[VizReadsArgs] {
  val companion: ADAMCommandCompanion = VizReads

  override def run(sc: SparkContext, job: Job): Unit = {
    val proj = Projection(contig, readMapped, readName, start, end)
    if (args.readPath.endsWith(".bam") || args.readPath.endsWith(".sam") || args.readPath.endsWith(".align.adam")) {
      VizReads.reads = sc.loadAlignments(args.readPath, projection = Some(proj))
      VizReads.readsRefName = args.readsRefName
    }

    if (args.variantsPath.endsWith(".vcf") || args.variantsPath.endsWith(".gt.adam")) {
      VizReads.variants = sc.loadGenotypes(args.variantsPath, projection = Some(proj))
      VizReads.variantsRefName = args.variantsRefName
    }

    val server = new org.eclipse.jetty.server.Server(args.port)
    val handlers = new org.eclipse.jetty.server.handler.ContextHandlerCollection()
    server.setHandler(handlers)
    handlers.addHandler(new org.eclipse.jetty.webapp.WebAppContext("mango-cli/src/main/webapp", "/"))
    server.start()
    println("View the visualization at: " + args.port)
    println("Frequency visualization at: /freq")
    println("Overlapping reads visualization at: /reads")
    println("Variant visualization at: /variants")
    println("Overall visualization at: /overall")
    server.join()
  }
}
