/******************************************************************************
 * KaHIPWrapper.java
 *
 * Example wrapper for Java integration of KaHIP via JNI
 *
 *****************************************************************************/

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.*;

public class KaHIPWrapper {

	private final static String NativeLibraryName = "wrapkahip";

	public final static int KAHIP_FAST         = 0;
	public final static int KAHIP_ECO          = 1;
	public final static int KAHIP_STRONG       = 2;
	public final static int KAHIP_FASTSOCIAL   = 3;
	public final static int KAHIP_ECOSOCIAL    = 4;
	public final static int KAHIP_STRONGSOCIAL = 5;

	static {
		try {
			System.loadLibrary(NativeLibraryName);
		} catch (UnsatisfiedLinkError e) {
			System.out.println("KaHIPWrapper - Error cannot find " + NativeLibraryName);
			e.printStackTrace();
		}
	}

	/*
	 * void kaffpa( int *n, int *vwgt, int *xadj, int *adjcwgt, int *adjncy, int
	 * *nparts, double *imbalance, bool suppress_output, int seed, int mode, int
	 * *edgecut, int *part);
	 */

	private static native void cnativeKaffpa(int n, int[] vwgt, int[] xadj,
			int[] adjcwgt, int[] adjncy, int nparts, double imbalance,
			boolean suppress_output, int seed, int mode,
			KaHIPWrapperResult returnInfo);

	public static KaHIPWrapperResult kaffpa(int n, int[] vwgt, int[] xadj,
			int[] adjcwgt, int[] adjncy, int nparts, double imbalance,
			boolean suppress_output, int seed, int mode) {

		KaHIPWrapperResult result = new KaHIPWrapperResult();

		cnativeKaffpa(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance,
				suppress_output, seed, mode, result);

		return result;
	}

	public static void main(String[] args) {
		// example: partition the graph in figure 4 of the KaHIP documentation
		int n                   = 5;
		int[] vwgt              = new int[0];
		int[] xadj              = new int[] {0, 2, 5, 7, 9, 12};
		int[] adjcwgt           = new int[0];
		int[] adjncy            = new int[] {1, 4, 0, 2, 4, 1, 3, 2, 4, 0, 1, 3};
		int nparts              = 2;
		double imbalance        = 0.4;
		boolean suppress_output = false;
		int seed                = 123456;
		int mode                = KaHIPWrapper.KAHIP_STRONG;

		KaHIPWrapperResult result = KaHIPWrapper.kaffpa(n, vwgt, xadj, adjcwgt, adjncy, 
								nparts, imbalance, suppress_output, seed, mode);

		System.out.println("=======\nResult:\n=======");

		System.out.println("edce cut = " + result.getEdgecut());
		System.out.println("Partitions:");
		for (int i=0; i<result.getPart().length; i++) 
			System.out.println("\tNode " + i + " belongs to block " + result.getPart()[i]);
	}

}
