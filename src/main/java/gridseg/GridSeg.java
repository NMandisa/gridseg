package gridseg;

import java.awt.Color;
import java.awt.Paint;
import java.awt.Point;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.fapesp.myutils.MyUtils;

public class GridSeg {

	private static final boolean DEBUG = false;
	private static final boolean GUIDEBUG = false;
	
	private static JFreeChart chart;
	private static ChartFrame frame;
	private static int keyCounter = 0;
	
	private static void plotBreaks(double[] breaksX, double[] breaksY, double minX, double maxX, double minY, double maxY) {	
		
		keyCounter = chart.getXYPlot().getDataset().getSeriesCount();
		XYSeriesCollection set = (XYSeriesCollection) chart.getXYPlot().getDataset();
		
		for (int i = 0; i < breaksX.length; i++) {			
			// plot x break:			
			XYSeries xy1 = new XYSeries(keyCounter);
			xy1.add(breaksX[i], minY);
			xy1.add(breaksX[i], maxY);			
			
			set.addSeries(xy1);
			chart.getXYPlot().getRenderer().setSeriesPaint(keyCounter, Color.black);
			chart.getXYPlot().getRenderer().setSeriesVisibleInLegend(keyCounter, false);
			keyCounter++;
			
			// plot y break:
			XYSeries xy2 = new XYSeries(keyCounter);
			xy2.add(minX, breaksY[i]);
			xy2.add(maxX, breaksY[i]);			
			
			set.addSeries(xy2);
			chart.getXYPlot().getRenderer().setSeriesPaint(keyCounter, Color.black);
			chart.getXYPlot().getRenderer().setSeriesVisibleInLegend(keyCounter, false);
			keyCounter++;
		}
		
		frame.pack();
		frame.repaint();
	}

	/**
	 * @param pts
	 * @param labels
	 */
	private static void plotMain(double[][] pts, int[] labels) {
		
		List<Paint> colors = Arrays.asList(ChartColor.createDefaultPaintArray());
		Collections.reverse(colors);

		XYSeriesCollection dataset = new XYSeriesCollection();
		
		HashSet<Integer> clusters = MyUtils.getUniqueElements(labels);
		Object[] uclusters = clusters.toArray();
		
		for (int i = 0; i < uclusters.length; i++) {
			int mi = (Integer) uclusters[i];
			String label = mi != 0 ? "Cluster" + mi : "Noise";
			
			XYSeries p = new XYSeries(label);
			for (int j = 0; j < labels.length; j++)
				if (labels[j] == mi)
					p.add(pts[j][0], pts[j][1]);
			
			dataset.addSeries(p);
		}
		
		// create chart:
		chart = ChartFactory.createXYLineChart("Plot", "x", "y", dataset);
		
		XYPlot xypl = chart.getXYPlot();
		xypl.setDomainGridlinesVisible(false);
		xypl.setRangeGridlinesVisible(false);
		NumberAxis na = (NumberAxis) xypl.getRangeAxis();
		na.setAutoRangeIncludesZero(false);
		
		XYLineAndShapeRenderer xr = (XYLineAndShapeRenderer) chart.getXYPlot().getRenderer();
		
		Shape circle = new Ellipse2D.Double(0, 0, 3, 3);
		
		for (int i = 0; i < dataset.getSeriesCount(); i++) {			
			xypl.getRenderer().setSeriesPaint(i, colors.get(i % colors.size()));
			xr.setSeriesShape(i, circle);
			xr.setSeriesLinesVisible(i, false);
			xr.setSeriesShapesVisible(i, true);
		}
		
//		chart.removeLegend();

		frame = new ChartFrame("Plot window", chart);
		frame.pack();
		frame.setDefaultCloseOperation(ChartFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}

	/**
	 * Grid segmentation and clustering algorithm. A cell is considered relevant
	 * if it has at least threshold points. Adjoining relevant cells are merged
	 * together with the same label.
	 * 
	 * @param data
	 *            examples, one per row
	 * @param nbreaks
	 *            number of grid divisions per dimension
	 *            
	 * @param minpts
	 * 			minimum number of points for a cluster to be regarded as relevant
	 * 
	 * @return cluster labels
	 */
	public static int[] gridseg(double[][] data, int nbreaks, int minpts) {

		int N = data.length;
		int m = data[0].length;

		int[] labels = new int[N];

		if (m != 2)
			throw new RuntimeException("Impl only for 2d data!");
		
		int threshold = Integer.MAX_VALUE;

		RealMatrix rm = new Array2DRowRealMatrix(data);

		double minX, maxX, minY, maxY;

		DescriptiveStatistics ds0 = new DescriptiveStatistics(rm.getColumn(0));
		DescriptiveStatistics ds1 = new DescriptiveStatistics(rm.getColumn(1));

		minX = ds0.getMin();
		maxX = ds0.getMax();
		minY = ds1.getMin();
		maxY = ds1.getMax();

		double[] breaksX = MyUtils.computeBreaks(minX, maxX, nbreaks);
		double[] breaksY = MyUtils.computeBreaks(minY, maxY, nbreaks);

		HashMap<Point, List<Integer>> gridmap = new HashMap<Point, List<Integer>>();
		for (int i = 0; i <= nbreaks; i++)
			for (int j = 0; j <= nbreaks; j++)
				gridmap.put(new Point(i, j), new ArrayList<Integer>());

		// ***************************************************
		// 1st step - get freq counts for each grid square
		// ***************************************************

		int[][] counts = new int[nbreaks + 1][nbreaks + 1];
		int[][] lbls = new int[nbreaks + 1][nbreaks + 1]; // labels
		boolean[][] visited = new boolean[nbreaks + 1][nbreaks + 1]; 

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

//		if (DEBUG) {
//			System.out.println("Counts matrix:");
//			MyUtils.print_matrix(counts);
//		}
		
		// determine the mean number of points in all cells:
		DescriptiveStatistics cellstats = new DescriptiveStatistics();
		for (int i = 0; i <= nbreaks; i++)
			for (int j = 0; j <= nbreaks; j++)
				cellstats.addValue(counts[i][j]);
		
		// set threshold:
		threshold = (int) Math.ceil(cellstats.getMean() + cellstats.getStandardDeviation());
		
		if (DEBUG)
			System.out.println("Using threshold = " + threshold);
		

		// ***************************************************
		// 2nd step - place labels and merge neighbors
		// ***************************************************
		int K = 0;
		
		for (int i = 0; i <= nbreaks; i++) {
			for (int j = 0; j <= nbreaks; j++) {
				if (visited[i][j])
					continue;

				if (counts[i][j] < threshold) {
					visited[i][j] = true;
					continue;
				}
				
				lbls[i][j] = ++K;
				
				// get relevant neighbors for grid cell (i,j):
				ArrayDeque<Point> neigh = getRelevantNeighborsList(i, j, counts, visited, threshold);
				while (neigh.size() > 0) {
					p = neigh.removeFirst();
					if (visited[p.x][p.y]) continue;
					
					lbls[p.x][p.y] = K;
					visited[p.x][p.y] = true;
					neigh.addAll(getRelevantNeighborsList(p.x, p.y, counts, visited, threshold));
				}
			}
		}
		
		// ***************************************************
		// 3rd step - remove outliers
		// ***************************************************

		// we consider as outliers cells that are relevant but have no
		// other adjoining relevant cells.
		
		for (int i = 0; i <= nbreaks; i++) {
			for (int j = 0; j <= nbreaks; j++) {
				if(!hasRelevantNeighbor(i, j, counts, threshold))
					lbls[i][j] = 0;
			}
		}
		
		// now actually assign the labels obtained:
		for (int i = 0; i <= nbreaks; i++)
			for (int j = 0; j <= nbreaks; j++) {
				// get the label for this cell
				int _lbl = lbls[i][j];
				p.x = i; p.y = j;
				// get the points for this cell:
				List<Integer> pts = gridmap.get(p);
				for (int pt : pts)
					labels[pt] = _lbl;
			}
		
		// we also consider as outliers clusters that have less than
		// minpts
		HashMap<Integer, Integer> ccounts = MyUtils.computeNumberOfOccurrences(labels);
		
		DescriptiveStatistics clusterstats = new DescriptiveStatistics();
		for (int val : ccounts.values())
			clusterstats.addValue(val);
		
		if (DEBUG)
			System.out.println("Using minpts = " + minpts);
		
		for (int i = 0; i < labels.length; i++) {
			int l = labels[i];
			if (l == 0) continue;
			if (ccounts.get(l) < minpts)
				labels[i] = 0;
		}
				
		
		// normalize labels as follows:
		// 0 - noise
		// 1 - first cluster
		// 2 - second cluster, etc.
		
		HashSet<Integer> uclusters = MyUtils.getUniqueElements(labels);
		int[] newlbls = new int[labels.length];
		Object[] ucl = uclusters.toArray();
		K = 0;
		for (int i = 0 ; i < ucl.length; i++) {
			int cl = (int) ucl[i];
			if (cl == 0) continue;
			K++;
			for (int j = 0; j < newlbls.length; j++)
				if(labels[j] == cl)
					newlbls[j] = K;
		}
		labels = newlbls;
		
		if (GUIDEBUG) {
			plotMain(data, labels);
			plotBreaks(breaksX, breaksY, minX, maxX, minY, maxY);
		}
		
		return labels;
	}
	
	private static boolean hasRelevantNeighbor(int x, int y, int[][] counts, double threshold) {
		
		int N = counts.length;

		final int[] offs = new int[] { -1, 0, 1 }; // offsets
		int nx, ny; // neighbor x and y coordinates

		for (int xo : offs) {
			for (int yo : offs) {

				if (xo == 0 && yo == 0)
					continue;

				nx = x + xo;
				ny = y + yo;

				if (nx < 0 || ny < 0 || nx >= N || ny >= N)
					continue;
				
				if (counts[nx][ny] >= threshold)
					return true;
			}
		}

		return false;
	}	

	private static ArrayDeque<Point> getRelevantNeighborsList(int x, int y, int[][] counts, boolean[][] visited, double threshold) {
		ArrayDeque<Point> neigh = new ArrayDeque<Point>();
		
		int N = counts.length;

		final int[] offs = new int[] { -1, 0, 1 }; // offsets
		int nx, ny; // neighbor x and y coordinates

		for (int xo : offs) {
			for (int yo : offs) {

				if (xo == 0 && yo == 0)
					continue;

				nx = x + xo;
				ny = y + yo;

				if (nx < 0 || ny < 0 || nx >= N || ny >= N)
					continue;
				
				if (!visited[nx][ny] && counts[nx][ny] >= threshold)
					neigh.add(new Point(nx, ny));
			}
		}

		return neigh;
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

		 File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/Y2.csv"); // topology data set
//		File csvData = new File("C:/Users/Cássio/workspace/gridseg/Y2cut.csv"); // topology data set
//		 File csvData = new	 File("C:/Users/Cássio/workspace/darbseg/resources/sample.csv"); //	 example given in the paper
//		 File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/sample3.csv"); //	 gaussians
//		 File csvData = new File("C:/Users/Cássio/workspace/darbseg/resources/simple.csv");
		 //		File csvData = new File("C:/Users/Cássio/workspace/gridseg/gridtest.csv");

		CSVParser parser = CSVParser.parse(csvData, Charset.defaultCharset(),
				CSVFormat.RFC4180.withDelimiter(' ').withSkipHeaderRecord().withIgnoreEmptyLines());

		List<CSVRecord> list = parser.getRecords();

		double[][] data = new double[list.size()][list.get(0).size()];

		for (int i = 0; i < list.size(); i++) {
			CSVRecord r = list.get(i);
			for (int j = 0; j < r.size(); j++)
				data[i][j] = Double.parseDouble(r.get(j));
		}

		int[] labels = gridseg(data, 35, 50);
		System.out.println("final partition:");
		MyUtils.print_array(labels);
	}

}
