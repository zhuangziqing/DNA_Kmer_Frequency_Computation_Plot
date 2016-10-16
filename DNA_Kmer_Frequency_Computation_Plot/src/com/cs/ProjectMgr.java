package com.cs;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;

public class ProjectMgr {
	/**
	 * @param args
	 */
	//Define four hashmap to store kmer Frequency when k=1,3,5,7
	static Map<String, Double> kmerFrequency1 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequency3 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequency5 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequency7 = new HashMap<String, Double>();
	//Define four hashmap to store kmer complement Frequency when k=1,3,5,7
	static Map<String, Double> kmerRevFrequency1 = new HashMap<String, Double>();
	static Map<String, Double> kmerRevFrequency3 = new HashMap<String, Double>();
	static Map<String, Double> kmerRevFrequency5 = new HashMap<String, Double>();
	static Map<String, Double> kmerRevFrequency7 = new HashMap<String, Double>();
	//define four hashmap to store kmer when kmer slide k bases instead of 1
	static Map<String, Double> kmerFrequencySlideK1 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequencySlideK3 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequencySlideK5 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequencySlideK7 = new HashMap<String, Double>();
	//define four hashmap to store kmer's reverse complement when kmer slide k bases instead of 1
	static Map<String, Double> kmerFrequencyRevSlideK1 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequencyRevSlideK3 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequencyRevSlideK5 = new HashMap<String, Double>();
	static Map<String, Double> kmerFrequencyRevSlideK7 = new HashMap<String, Double>();
	
	public static void main(String[] args) {
		int k1 = 1;
		int k3 = 3;
		int k5 = 5;
		int k7 = 7;
		
		String filePath = "hs_ref_GRCh38_chr17.fa";
		ProjectMgr mgr = new ProjectMgr();
		/*kmerFrequency1 = bioinformatcisHW3.getFreq(filePath, k1);
		kmerRevFrequency1 = bioinformatcisHW3.getRevComp(kmerFrequency1);
		BioinformaticsHW3.getImage(kmerRevFrequency1);
		*/
		kmerFrequencySlideK7 = mgr.getFreqKSlides(filePath, k7);
		kmerFrequencyRevSlideK7 = mgr.getRevComp(kmerFrequencySlideK7);
		ProjectMgr.getImage(kmerFrequencyRevSlideK7);
	}

	public Map<String, Double> getFreq(String filePath, int k){
		
		BufferedReader br = null;
		FileReader fr = null;
		String currentLine = "";
		Map<String, Integer> kmerCount = new HashMap<String, Integer>();
		Map<String, Double> kmerFreq = new HashMap<String, Double>();
		
		try {
			fr = new FileReader(filePath);
			br = new BufferedReader(new FileReader(filePath));
			String s = br.readLine();
			while(s != null){
				if(s.contains("gi")){
					s = br.readLine();
					continue;
				}else{
					currentLine = currentLine + s;
					int currentLineLength = currentLine.length();
					int positionLine = 0;
					
					while(positionLine <= currentLineLength - k){
						String kmer = currentLine.substring(positionLine, positionLine + k);
						
						if(kmerCount.containsKey(kmer)){
							kmerCount.put(kmer, kmerCount.get(kmer) + 1);
						}else{
							kmerCount.put(kmer, 1);
						}
						positionLine++;
					}
					currentLine = currentLine.substring(currentLineLength - k);
					s = br.readLine();
				}
			}
		}catch(FileNotFoundException e){
			e.printStackTrace();
		}catch(IOException e){
			e.printStackTrace();
		}finally{
			try {
				fr.close();
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		Set<Map.Entry<String, Integer>> set = kmerCount.entrySet();
		int totalNum = 0;
		//use this for loop to get the total number of all kmer
		for(Iterator<Map.Entry<String, Integer>> it = set.iterator(); it.hasNext(); ){
			Map.Entry<String, Integer> mapEntry = it.next();
			totalNum = totalNum + mapEntry.getValue();
			//System.out.println("When Kmer = " + mapEntry.getKey() + ", its nubmer = " + mapEntry.getValue());
		}
		//System.out.println("total number " + totalNum);
		//use this for loop to go through all the elements in hashmap and caculate the frequency
		for(Iterator<Map.Entry<String, Integer>> it = set.iterator(); it.hasNext(); ){
			Map.Entry<String, Integer> mapEntry = it.next();
			double value = (double)mapEntry.getValue()/(double)totalNum;
			kmerFreq.put(mapEntry.getKey(), value);
			//System.out.println("key = " + mapEntry.getKey() + ", value = " + value);
		}
		return kmerFreq;
	}

	//this function input original mapFrequency and output reverse complement frequency
	public Map getRevComp(Map<String, Double> mapOriginal){
		Map<String, Double> mapRevComp = new HashMap<String, Double>();
		Set<Map.Entry<String, Double>> set = mapOriginal.entrySet();
		for(Iterator<Map.Entry<String, Double>> it = set.iterator(); it.hasNext(); ){
			Map.Entry<String, Double> mapEntry = it.next();
			StringBuffer temp = new StringBuffer(mapEntry.getKey());
			String revComp;
			for(int i=0; i<temp.length(); i++){
				
				if(String.valueOf(temp.charAt(i)).equals("A")){
					temp.setCharAt(i, 'T');
				}else if(String.valueOf(temp.charAt(i)).equals("T")){
					temp.setCharAt(i, 'A');
				}else if(String.valueOf(temp.charAt(i)).equals("C")){
					temp.setCharAt(i, 'G');
				}else if(String.valueOf(temp.charAt(i)).equals("G")){
					temp.setCharAt(i, 'C');
				}
			}
			revComp = temp.reverse().toString();
			mapRevComp.put(revComp, mapEntry.getValue());
		}
		
		Set<Map.Entry<String, Double>> setFrequency = mapRevComp.entrySet();
		for(Iterator<Map.Entry<String, Double>> it = setFrequency.iterator(); it.hasNext(); ){
			Map.Entry<String, Double> mapEntry = it.next();
			//System.out.println("kmer's reverse complement is " + mapEntry.getKey() + ", frequency is " + mapEntry.getValue() + "  ");
		}
		return mapRevComp;
	}
	
	//this function calcuate frequency when kmer slide every k bases
	public Map<String, Double> getFreqKSlides(String filePath, int k){
		
		BufferedReader br = null;
		FileReader fr = null;
		String currentLine = "";
		Map<String, Integer> kmerCount = new HashMap<String, Integer>();
		Map<String, Double> kmerFreqKSlides = new HashMap<String, Double>();
		
		try {
			fr = new FileReader(filePath);
			br = new BufferedReader(new FileReader(filePath));
			String s = br.readLine();
			while(s != null){
				if(s.contains("gi")){
					s = br.readLine();
					continue;
				}else{
					currentLine = currentLine + s;
					int currentLineLength = currentLine.length();
					int positionLine = 0;
					
					while(currentLineLength - positionLine >= k){
						String kmer = currentLine.substring(positionLine, positionLine + k);
						
						if(kmerCount.containsKey(kmer)){
							kmerCount.put(kmer, kmerCount.get(kmer) + 1);
						}else{
							kmerCount.put(kmer, 1);
						}
						positionLine = positionLine + k;
					}
					currentLine = currentLine.substring(currentLineLength - k);
					s = br.readLine();
				}
			}
		}catch(FileNotFoundException e){
			e.printStackTrace();
		}catch(IOException e){
			e.printStackTrace();
		}finally{
			try {
				fr.close();
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		Set<Map.Entry<String, Integer>> set = kmerCount.entrySet();
		int totalNum = 0;
		//use this for loop to get the total number of all kmer
		for(Iterator<Map.Entry<String, Integer>> it = set.iterator(); it.hasNext(); ){
			Map.Entry<String, Integer> mapEntry = it.next();
			totalNum = totalNum + mapEntry.getValue();
			//System.out.println("When Kmer = " + mapEntry.getKey() + ", its nubmer = " + mapEntry.getValue());
		}
		//System.out.println("total number " + totalNum);
		//use this for loop to go through all the elements in hashmap and caculate the frequency
		for(Iterator<Map.Entry<String, Integer>> it = set.iterator(); it.hasNext(); ){
			Map.Entry<String, Integer> mapEntry = it.next();
			double value = (double)mapEntry.getValue()/(double)totalNum;
			kmerFreqKSlides.put(mapEntry.getKey(), value);
			//System.out.println("key = " + mapEntry.getKey() + ", value = " + value);
		}
		return kmerFreqKSlides;
	}
	
	public static void getImage(Map kmerRevFrequency){
		CategoryDataset dataset = getDataSet(kmerRevFrequency);
		JFreeChart chart = ChartFactory.createBarChart3D(
							"K-mer Reverse Complement Frequency (Chr 17)", // histogram title
							"K-mer (k = 1)", // x-axis
							"Frequency", // y-axis
							dataset, // dataset
							PlotOrientation.VERTICAL, // direction: vertical or horizontal
							true, 	// show the image or not(for simple histogram must be false)
							false, 	// produce tool
							false 	// produce url or not
							);
							
		FileOutputStream kmerImage = null;
		try {
			kmerImage = new FileOutputStream("D://kmer.jpg");
			ChartUtilities.writeChartAsJPEG(kmerImage,0.5f,chart,1600,1200,null);
		} catch(IOException ex){
			ex.printStackTrace();
		} catch(Exception exc){
			exc.printStackTrace();
		} finally {
			try {
				if(kmerImage != null){
				kmerImage.close();
				}
			} catch (Exception e) {e.printStackTrace();}
		}
	}
	
	private static CategoryDataset getDataSet(Map k) {
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		//dataset.addValue(100, "", "apple");
		Set<Map.Entry<String, Double>> set = k.entrySet();
		for(Iterator<Map.Entry<String, Double>> it = set.iterator(); it.hasNext(); ){
			Map.Entry<String, Double> mapEntry = it.next();
			dataset.addValue(mapEntry.getValue(), "", mapEntry.getKey());
			System.out.println("key = " + mapEntry.getKey() + ", value = " + mapEntry.getValue());
		}
		return dataset;
	}
	
}
