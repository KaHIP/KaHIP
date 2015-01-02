/******************************************************************************
 * KaHIPWrapper.cpp
 *
 * Example wrapper for Java integration of KaHIP via JNI
 *
 ******************************************************************************
 * Copyright (C) 2014 Uniserv GmbH
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/


#include "KaHIPWrapper.h"
#include "kaHIP_interface.h"
#include <stdlib.h>

/*
	 void kaffpa( n, int *vwgt, int *xadj, int *adjcwgt, int *adjncy, int
	 nparts, double *imbalance, bool suppress_output, int seed, int mode, int
	 *edgecut, int *part);
*/

JNIEXPORT void JNICALL Java_KaHIPWrapper_cnativeKaffpa
  (JNIEnv *pEnv, jclass jObj, jint jn, jintArray jvwgt, jintArray jxadj, jintArray jadjcwgt, jintArray jadjncy,
  jint jnparts, jdouble jimbalance, jboolean jsuppress_output, jint jseed, jint jmode, jobject jRetObj)
{
	jsize i;
	jsize len;
	jfieldID fid;
	jclass resultClass;
	jintArray jpart;
	jint jedgecut;
	jint *vwgtJArray, *xadjJArray, *adjcwgtJArray, *adjncyJArray, *partJArray;
	int n, nparts, seed, mode;
	int *vwgt, *xadj, *adjcwgt, *adjncy, *part;
	double imbalance;
	bool suppress_output;
	int edgecut;

	// convert java input data types to C++ data types
	n = (jint)jn;

	len = pEnv->GetArrayLength(jvwgt);
	if (len == 0) {
		vwgt = 0;
	} else {
		vwgtJArray = pEnv->GetIntArrayElements(jvwgt, 0);
		vwgt = (int*)malloc(len * sizeof(int));
		for (i = 0; i < len; i++) vwgt[i] = (int)vwgtJArray[i];
	}

	len        = pEnv->GetArrayLength(jxadj);
	xadjJArray = pEnv->GetIntArrayElements(jxadj, 0);
	xadj       = (int*)malloc(len * sizeof(int));

	for (i = 0; i < len; i++) {
                xadj[i] = (int)xadjJArray[i];
        }

	len = pEnv->GetArrayLength(jadjcwgt);
	if (len == 0) {
		adjcwgt = 0;
	} else {
		adjcwgtJArray = pEnv->GetIntArrayElements(jadjcwgt, 0);
		adjcwgt = (int*)malloc(len * sizeof(int));
		for (i = 0; i < len; i++) adjcwgt[i] = (int)adjcwgtJArray[i];
	}

	len          = pEnv->GetArrayLength(jadjncy);
	adjncyJArray = pEnv->GetIntArrayElements(jadjncy, 0);
	adjncy       = (int*)malloc(len * sizeof(int));

	for (i = 0; i < len; i++) {
                adjncy[i] = (int)adjncyJArray[i];
        }

	nparts          = (int)jnparts;
	imbalance       = (double)jimbalance;
	suppress_output = (bool)jsuppress_output;
	seed            = (int)jseed;
	mode            = (int)jmode;

	// allocate memory for output parameters of kaffpa
	part = (int*)malloc(n * sizeof(int));

	// KaHIP call for KaFFPa graph partitioner
	kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, suppress_output, seed, mode, &edgecut, part);

	// convert output to JNI objects
	partJArray = (jint *)malloc(n * sizeof(jint));

	for (int i=0; i<n; i++) {
		partJArray[i] = (jint)part[i];
        }

	jpart = pEnv->NewIntArray(n);
	pEnv->SetIntArrayRegion(jpart, 0, n, (const jint *)partJArray);

	jedgecut = (jint)edgecut;

	// populate KaHIPWrapperResult object
	resultClass = pEnv->GetObjectClass(jRetObj);

	fid = pEnv->GetFieldID(resultClass, "edgecut", "I");
	pEnv->SetIntField(jRetObj, fid, jedgecut);

	fid = pEnv->GetFieldID(resultClass, "part", "[I");
	pEnv->SetObjectField(jRetObj, fid, jpart);

	// free memory
	free(vwgt); pEnv->ReleaseIntArrayElements(jvwgt, vwgtJArray, 0);
	free(xadj); pEnv->ReleaseIntArrayElements(jxadj, xadjJArray, 0);
	free(adjcwgt); pEnv->ReleaseIntArrayElements(jadjcwgt, adjcwgtJArray, 0);
	free(adjncy); pEnv->ReleaseIntArrayElements(jadjncy, adjncyJArray, 0);
	free(part); pEnv->ReleaseIntArrayElements(jpart, partJArray, 0);
}
