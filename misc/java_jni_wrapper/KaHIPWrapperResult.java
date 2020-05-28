/******************************************************************************
 * KaHIPWrapperResult.java
 *
 * Example wrapper for Java integration of KaHIP via JNI
 *****************************************************************************/


public class KaHIPWrapperResult {
	private int edgecut;
	private int[] part;
	
	public int getEdgecut() {
		return edgecut;
	}
	public void setEdgecut(int edgecut) {
		this.edgecut = edgecut;
	}
	public int[] getPart() {
		return part;
	}
	public void setPart(int[] part) {
		this.part = part;
	}
	
}
