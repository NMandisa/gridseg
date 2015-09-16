package gridseg;

import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import br.fapesp.myutils.MyUtils;

public class GridSeg {
	
	private static final boolean DEBUG = true;
	
	public static int[] gridseg(double[][] data, int nbreaks, int threshold) {
		
		int N = data.length;
		int m = data[0].length;
		
		int[] labels = new int[N];
		
		if (m != 2)
			throw new RuntimeException("Impl only for 2d data!");
		
		RealMatrix rm = new Array2DRowRealMatrix(data);
		
		double minX, maxX, minY, maxY;
		
		DescriptiveStatistics ds0 = new DescriptiveStatistics(rm.getColumn(0));
		DescriptiveStatistics ds1 = new DescriptiveStatistics(rm.getColumn(1));
		
		minX = ds0.getMin(); maxX = ds0.getMax();
		minY = ds1.getMin(); maxY = ds1.getMax();
		
		double[] breaksX = MyUtils.computeBreaks(minX, maxX, nbreaks);
		double[] breaksY = MyUtils.computeBreaks(minY, maxY, nbreaks);
		
		HashMap<Point, List<Integer>> gridmap = new HashMap<Point, List<Integer>>();
		for (int i = 0; i <= nbreaks; i++)
			for (int j = 0; j <= nbreaks; j++)
				gridmap.put(new Point(i,j), new ArrayList<Integer>());
		
		// ***************************************************
		// 1st step - get freq counts for each grid square 
		// ***************************************************
		
		int[][] counts = new int[nbreaks+1][nbreaks+1];
		int[][] lbls = new int[nbreaks+1][nbreaks+1]; // labels
		
		Point p = new Point();
		
		for (int i = 0; i < N; i++) {
			// find the grid position the i-th point falls on
			int gridX = getGridPos(data[i][0], breaksX, true);
			int gridY = getGridPos(data[i][1], breaksY, false);
			p.x = gridY; // because y refers to the rows
			p.y = gridX; // because x refers to the columns
			List<Integer> pts = gridmap.get(p);
			pts.add(i);			
			counts[gridY][gridX]++;
		}		
		

		if (DEBUG) {
			System.out.println("Counts matrix:");
			MyUtils.print_matrix(counts);
		}		
		
		// ***************************************************
		// 2nd step - place labels and merge neighbors
		// ***************************************************		
		int K = 1;
		
		for (int i = 0; i <= nbreaks; i++)
			for (int j = 0; j <= nbreaks; j++)
				if (counts[i][j] >= threshold) {
					// 
				}
					
		
		
		
		return labels;		
	}
	
	
	  private static int getGridPos(double p, double[] breaks, boolean isXdim) {
		  int pos = -1;
		  
		  for (int i = 0; i < breaks.length; i++)
			  if (p <= breaks[i]) {
				  pos = i;
				  break;
			  }
		  
		  if (pos == -1) // point falls after the last break
			  pos = breaks.length;
		  
		  if (!isXdim)
			  pos = breaks.length - pos;
		  
		  return pos;
	}


	public static void main(String args[]) throws IOException {
	    	
//			File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/Y2.csv"); // topology data set
//	    	File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/sample.csv"); // example given in the paper
//	    	File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/sample3.csv"); // gaussians
//	    	File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/simple.csv");
			File csvData = new File("C:/Users/Cássio/workspace/gridseg/gridtest.csv");
	    	
	        CSVParser parser = CSVParser.parse(csvData, Charset.defaultCharset(), CSVFormat.RFC4180.withDelimiter(' ').withSkipHeaderRecord().withIgnoreEmptyLines());
	        
	        List<CSVRecord> list = parser.getRecords();
	        
	        double[][] data = new double[list.size()][list.get(0).size()];
	        
	        for (int i = 0; i < list.size(); i++) {
	        	CSVRecord r = list.get(i);
	        	for (int j = 0; j < r.size(); j++)
	        		data[i][j] = Double.parseDouble(r.get(j));
	        }
	        
	      gridseg(data, 2, 3);
	        
	    }
	    

}
