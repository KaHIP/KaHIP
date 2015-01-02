/******************************************************************************
 * KaHIPWrapperResult.java
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
