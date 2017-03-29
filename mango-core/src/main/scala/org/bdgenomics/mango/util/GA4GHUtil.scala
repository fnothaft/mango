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
package org.bdgenomics.mango.core.util

import ga4gh.Variants.Variant

object GA4GHUtil {
  def main(args: Array[String]): Unit = {
    println("Hello Justin")
    val x_builder = ga4gh.Variants.Variant.newBuilder()
    x_builder.setEnd(342424)
    val x = x_builder.build()

    val myString = com.google.protobuf.util.JsonFormat.printer().print(x)
    println("myString: " + myString)

  }
}

