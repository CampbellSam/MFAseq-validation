import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.IllegalFormatException;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceRecord;

public class Compare {
	
	 SAMFileReader BSF_eS_reader; 	//filereader variables to hold SAMFileReaders passed in constructor
	 SAMFileReader BSF_G2_reader; 
	 SAMFileReader PCF_eS_reader; 
	 SAMFileReader PCF_G2_reader; 
	 
	 int BSF_fileOneTC = 0;			//variable to hold the total aligned read count for each file
	 int BSF_fileTwoTC = 0;		
	 int PCF_fileOneTC = 0;
	 int PCF_fileTwoTC = 0;
	 
	 BAMFile fileOneBSF = null;		//create BAMFile object for each file - allows counting methods to be performed
	 BAMFile fileTwoBSF = null;
	 BAMFile fileOnePCF = null;
	 BAMFile fileTwoPCF = null;
	 
	 Map <String, Integer> map = null;		//map to contain chromosomes and respective lengths
	 Map <String, Double> weights = null;   //map to contain chromosome  names and weights
	 ArrayList al = null;
	 int totalLength = 0;
	 
	
	 /*
	  * Single constructor for Compare class
	  * assigns the file readers passed to constructor - **error handling if file is null**
	  * gets total read count and adds chromosomes to map
	  */
	public Compare (SAMFileReader file1_phase1, SAMFileReader file1_phase2, SAMFileReader file2_phase1, SAMFileReader file2_phase2)
	{	
		BSF_eS_reader = file1_phase1; 
		BSF_G2_reader = file1_phase2;
		PCF_eS_reader = file2_phase1;
		PCF_G2_reader = file2_phase2;
		
		fileOneBSF = new BAMFile(BSF_eS_reader);
		BSF_fileOneTC = fileOneBSF.getTotalReadCount();		//get total read count for each file
		fileTwoBSF = new BAMFile(BSF_G2_reader);
		BSF_fileTwoTC = fileTwoBSF.getTotalReadCount();
		fileOnePCF = new BAMFile(PCF_eS_reader);
		PCF_fileOneTC = fileOnePCF.getTotalReadCount();
		fileTwoPCF = new BAMFile(PCF_G2_reader);
		PCF_fileTwoTC = fileTwoPCF.getTotalReadCount();
		
		map = new HashMap<String, Integer>();		//create map which will contain each chromosome and the corresponding length
		

		SAMFileHeader fh;	
		int numChrRecs = 0;			// int to contain number or chromosomes
		totalLength = 0;
		try
		{
			fh = BSF_eS_reader.getFileHeader();		//instantiate fileHeader object with first file SAMFileReader
			for (SAMSequenceRecord seq : fh.getSequenceDictionary().getSequences()) 
		    {
				map.put(seq.getSequenceName(),seq.getSequenceLength());		// add 'name' - chromosome and length to map
			    numChrRecs ++;				// increment number of chromosomes
			    totalLength += seq.getSequenceLength();
			    //System.out.println(seq.getSequenceName());
	   	    }		
		 }
		 catch(IllegalFormatException ife)
		 {
			 System.err.println("BAM/SAM incorrect file header format");
			 System.exit(0);
		 }
		
		 
		weights = new HashMap <String, Double>();		//new hash to contain chromosome name and weights
	   
		//for each chromosome - weight is %
		for (String key : map.keySet())			//loop through each chromosome and calculate weight
		{
			 int length = map.get(key);			//chromosome length
			 double weight = (double) length / totalLength * 100;		//weight = length / totalLength 
			 weights.put(key, weight);			//add name and weight to hash
			 //writer.write(key + "\t" + length + "\t" + weight + "\n");				 
		}
		
		al = new ArrayList();		//arraylist to contain chromosome names
		int constant = 10;
		for (String key : weights.keySet())		//for each chromosome in weights hash
		{
			for (int j = 0; j < weights.get(key) * constant; j++)		//add chr name the number of times equal to the weight * the constant
			{
				al.add(key);		//the number of times added is dependent on weight
			}
		}
	}
	
	/*
	 * method to return the length of a specified chromosome
	 * will return 0 if string is empty
	 */
	public int getChromosomeLength (String chromosome)
	{
		int length = 0;
		if (chromosome != null && chromosome != "")		//if variable is not empty or = 0 
		{	
			try
			{
				length = map.get(chromosome);		//get the length from the map and return this
				return length ;
			
			}
			catch(NullPointerException n)
			{
				n.printStackTrace();
				System.out.println("Could not find length - chromosome not found");
			}
		}
		else
		{
			return length;		//returns 0
		}
		return length;			//returns 0
	}
	
	
	/*
	 * Method to remove a chromosome from the map
	 * void - does not have a return value
	 */
	public void removeChromosome(String chromosome)
	{
		if (chromosome != null && chromosome != "")
		{
			map.remove(chromosome);
			al.removeAll(Collections.singleton(chromosome));
			//System.out.println("Removed: " + chromosome);
		}
		else
		{
			System.err.println("Could not remove comparison region from chromosome selection");
		}
	}
	
	
	/*
	 * Method to return a randomly selected chromosome
	 * Chromosome weights are determined to allow more accurate
	 * random sampling.
	 * Parameters - hash containing chromosomes and length
	 * + the total length of the genome
	 */
	public String randomChromosome()
	{
	
		/*
		//output file to contain weights
	    Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("Chromosome_weights.txt"), "utf-8"));
	    writer.write("Chromosome\tLength\tWeight	\n");
	    */
					
		Random rand = new Random();
	    int num = rand.nextInt(al.size());			//generate a new random int with the maximum being the array length
	    String chromosome = (String) al.get(num);	//return chromosome name present at random index of array		
	    System.out.println(chromosome);
	    /*
	    System.out.println("Whole array in randomChromosome(): ");
	    for (int i= 0; i < al.size(); i++)
	    {
	    	System.out.println(al.get(i));
	    }
	    */
	    
		return chromosome;
	}
	
	
	
	// need to add tests to ensure input is correct - no nulls or incorrect format etc.
	public double mannWhitneyTest(String chromosome, int start, int end, int window_size)
	{
		 double pvalue = 0.0;
		
		 ArrayList bsf = new ArrayList();		//array to hold list of ratios for each form
		 ArrayList pcf = new ArrayList();
		 
		 //start & end of 60kb region
		 int chr_start = start;		
		 int chr_end = end;
		 
		 //loop through this region, incrementing by window size
		 for (int j = 0; j < chr_end; j = j + window_size)
		 {
			 int window_end = chr_start + window_size;		//end of the window = start coordinate + window size
			 if (window_end - chr_start < 2499)				//if the window is smaller than the size wanted the loop breaks
			 {
				 break;
			 }
			 
			 int BSF_eS_R = fileOneBSF.getRegionCount(chromosome, chr_start, window_end);		//get the aligned read count for this segment for each file
			 int BSF_G2_R = fileTwoBSF.getRegionCount(chromosome, chr_start, window_end);
			 int PCF_eS_R = fileOnePCF.getRegionCount(chromosome, chr_start, window_end);
			 int PCF_G2_R = fileTwoPCF.getRegionCount(chromosome, chr_start, window_end);
				 
			 double BSFphase_ratio = fileOneBSF.getRatio(BSF_fileOneTC, BSF_eS_R, BSF_fileTwoTC, BSF_G2_R);		//call get Ratio method with numbers generated
			 double PCFphase_ratio = fileOnePCF.getRatio(PCF_fileOneTC, PCF_eS_R, PCF_fileTwoTC, PCF_G2_R);
			 
			 bsf.add(BSFphase_ratio);		//add current ratio to ArrayList
			 pcf.add(PCFphase_ratio);
				 
			 chr_start = window_end + 1;		//new start coordinate
			 if (chr_start > chr_end)
			 {
				 break;			//if the new start coordinate is out with the 60kb window
			 }
			 //System.out.println(BSFphase_ratio + "\t" + PCFphase_ratio);
		 }
		 double [] bsfComp_array = new double[bsf.size()];		//double array to hold ratio values
		 double [] pcfComp_array = new double[pcf.size()];
		 
		 for (int k = 0; k < bsfComp_array.length; k++) 		//cannot find a way to directly cast Double arrayList to double [] and avoid this loop
		 {
			 bsfComp_array[k] = (double) bsf.get(k);
			 pcfComp_array[k] = (double) pcf.get(k);
		 }
		
		 MannWhitneyUTest mwt = new MannWhitneyUTest();		//create mwt object - takes double [] as parameters
		 pvalue = mwt.mannWhitneyUTest(bsfComp_array, pcfComp_array);		//perform calculation using ratio arrays
		 //System.out.println("x: " + pvalue); 
		
		return pvalue;			//return p-value
	}

	
}
