import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.IllegalFormatException;
import java.util.Map;
import java.util.Random;


import net.sf.samtools.*;
import picard.*;

/*
 * Main class  - file IO
 * finds x
 * finds y (* iterations)
 * calculates result
 *
 * values for hardcoded variables 
 * may need to be changed before running
 */
public class calculateRatio {
	
	public static void main(String[] args) throws IOException {
		
		 File BSF_eS = null;			//File variables to contain file input
		 File BSF_eS_index = null;
		 File BSF_G2 = null;
		 File BSF_G2_index = null;
		
		 File PCF_eS = null;
		 File PCF_eS_index = null;
		 File PCF_G2 = null;
		 File PCF_G2_index = null;
			  
		 String randChr = "";			//random chromosome for finding y
		 int start = 0;					//start of 60kb region
		 int end = 0;					// end of region/chromosome
		 int iterations = 12000;			//user input?? number of iterations of y
		 int region = 59999;			//60kb region
		 int window_size = 2500;		//user input
		 window_size -= 1;
		 
		 String comp_chr = "Tb427_telo40_v2";		//current comparison region
		 
		 double x = 0;
		 double y = 0;
		 int counter = 0;		//increments if y >= x
		 double result = 0;		//contains result - counter / iterations
		 
		 //8 input files
		 BSF_eS = new File ("/Users/samanthacampbell/Desktop/bam/T.brucei427H9B89ADXX_BSF427-EarleyS_GTCCGC.bam"); 
		 BSF_eS_index = new File(BSF_eS.getPath() + ".bai");
		 BSF_G2 = new File ("/Users/samanthacampbell/Desktop/bam/T.brucei427H9B89ADXX_BSF427-G21M_TAGCTT.bam"); 
		 BSF_G2_index = new File(BSF_G2.getPath() + ".bai");
		 
		 PCF_eS = new File ("/Users/samanthacampbell/Desktop/bam/T.brucei427H9B89ADXX_PCF427-EarlyS_ACTTGA.bam");
		 PCF_eS_index = new File(PCF_eS.getPath() + ".bai");
		 PCF_G2 = new File ("/Users/samanthacampbell/Desktop/bam/T.brucei427H9B89ADXX_PCF427-G21M_GATCAG.bam");
		 PCF_G2_index = new File(PCF_G2.getPath() + ".bai");
		 
		 File file = new File ("TrypPCF_BSF_G2eS12000.bed");
		 Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		 
		 String track_def = "track name = \"PCF & BSF\"\tdescription = \"User Supplied Track\"";
		 //writer.write(track_def);
		 
		 // create SAMFileReader from each input file
		 SAMFileReader BSF_eS_reader = new SAMFileReader(BSF_eS, BSF_eS_index);			
		 SAMFileReader BSF_G2_reader = new SAMFileReader(BSF_G2, BSF_G2_index);
		 
		 SAMFileReader PCF_eS_reader = new SAMFileReader(PCF_eS, PCF_eS_index);
		 SAMFileReader PCF_G2_reader = new SAMFileReader(PCF_G2, PCF_G2_index);
		
		 // creates new Compare class object - performs 2.5kb iterations
		 // compare class creates BAMFile object in constructor
		 Compare compare = new Compare(BSF_eS_reader, BSF_G2_reader, PCF_eS_reader, PCF_G2_reader);		
		 
		 
		 /*
		  * Calculate x first, loop through in 2.5kb's - telo40 1->end
		  * remove last window if smaller than 2.5kb
		  */
		 start = 1;
		 end = compare.getChromosomeLength(comp_chr);		//end = length of chromosome
		 x = compare.mannWhitneyTest(comp_chr, start, end, window_size);	//pass start, end, chromosome and iteration window size - returns pvalue from mww
		 compare.removeChromosome(comp_chr);		//remove comparison chromosome from list so it cannot be chosen in random selection
		 
		 System.out.println("x: " + x);
		 
		 
		 //choose random 60kb region for the number of iterations
		 for(int i = 0; i < iterations; i++)
		 {
			 randChr = compare.randomChromosome();			//choose random chromosome
			 end = compare.getChromosomeLength(randChr);	//get chromosome length
			 		 	 
		 	 while (end < region)		//test if the length of the chromosome is >60kb, try again while it is not
		  	 {
		 		 randChr = compare.randomChromosome();
				 end = compare.getChromosomeLength(randChr);
		  	 }
		 	 
		 	 int start_max = end - region;				//maximum start of segment - last coordinate minus segment size
		 	 Random random = new Random();
		 	 start = random.nextInt(start_max);			//randomly generate start coordinate using maximum
			 int region_end = start + region;			// end of segment = start + size of segment
			 			 
			 y = compare.mannWhitneyTest(randChr, start, region_end, window_size);		//pass start, end, chromosome and iteration window size - returns pvalue from mww
			 //System.out.println("y: " + y);
			 			 
			 if (y <= x)			//test if y is less than or equal to x & increment counter if true
			 {
				 counter++;
				// writer.write("\n" + randChr + "\t" + start + "\t" + region_end + "\t" + y);
				 System.out.println(randChr + "\t" + start + "\t" + region_end + "\t" + y);
			 }
			 writer.write("\n" + randChr + "\t" + start + "\t" + region_end + "\t" + y);
		     System.out.println(i);
		 } 
		 
		 BigDecimal temp_count = new BigDecimal(counter);				// counter from integer to bigDecimal - allows calculation
		 BigDecimal temp_iter = new BigDecimal(iterations);				// same for number of iterations
		 BigDecimal temp_result = temp_count.divide(temp_iter, 6, RoundingMode.HALF_EVEN);		//	counter / iterations

		 result = temp_result.doubleValue();		//result as double
		 System.out.println("counter: " + counter + "\nresult: " + result);
		 writer.close();
		 		 
	}

}
