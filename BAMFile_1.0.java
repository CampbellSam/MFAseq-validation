import java.math.BigDecimal;
import java.math.RoundingMode;

import net.sf.samtools.*;
import picard.*;

public class BAMFile {
			
private SAMFileReader currentSam;

	//class constructor
	public BAMFile(SAMFileReader sam)		//SAMFileReader initialised with bam file and corresponding .bai index file
	{
		currentSam = sam;	
	}
	
	//method to return the number of aligned reads in the current SAM/BAM file
	public int getTotalReadCount()
	{
		
		AbstractBAMFileIndex abstractIndex = (AbstractBAMFileIndex) currentSam.getIndex();		//creates index object from SAMFileReader
		int nRefs = abstractIndex.getNumberOfReferences();										//number of references in index (chromosomes)
		
		int nReads = 0;						//int to contain the number of aligned reads
		for (int i = 0; i < nRefs; i++)
		{
			BAMIndexMetaData metaData = abstractIndex.getMetaData(i);		//create metadata object for current chromosome (i) from index
			int numRecs = metaData.getAlignedRecordCount();					//number of reads that have been aligned to the current chromosome, can also use getUnalignedRecordCount()
			nReads += numRecs;												// increment the total number of reads
		}
		
		return nReads;		//returns the total number of aligned reads
	}
	
	
	
	//method to return the total reads aligned to a region of a chromosome - the chr and region are passed as parameter
	public int getRegionCount(String chr, int start, int end)
	{				
		SAMFileReader temp = currentSam;
		SAMRecordIterator it = temp.queryContained(chr, start, end);		//create iterator object for region of chromosome contained within parameters
		
	    int counter = 0;			// counter will track the number of aligned reads
		while (it.hasNext())
		{
			SAMRecord sam = it.next();		//SAMRecord created with read from current iteration - MAY NOT BE NEEDED
			counter++;						// increment counter
		}	
		it.close();
		
		return counter;				// return number of aligned reads in region 
	}
	
	 
    /*
     * Calculate ratio of reads aligned to region between the two files 
     * (diff. cell cycle phases)
     *    regionCounts1 / regionCounts2 * totalCounts2 / totalCounts1
     * Each numerical value cast to big decimal for calculation   
     */
	public double getRatio(int tc1, int r1, int tc2, int r2)
	{
		double ratio = 0;
		BigDecimal temp_ratio = new BigDecimal(0);
		
		if (r1 == 0 || r2 == 0)		// if either file has zero coverage for the region, ratio = 0.001
	    {
			temp_ratio = new BigDecimal(0.001);
	    }
	    else
	    {   
		    BigDecimal total1 = new BigDecimal(tc1);	
		    BigDecimal total2 = new BigDecimal(tc2);
		     
		    BigDecimal region1 = new BigDecimal(r1);
		    BigDecimal region2 = new BigDecimal(r2);  
		     
		    BigDecimal top = region1.multiply(total2);
		    BigDecimal bottom = region2.multiply(total1);
		    temp_ratio = top.divide(bottom, 6, RoundingMode.HALF_EVEN);
		     
		  //System.out.println("top: " + top + "\tbottom: " + bottom + "\nratio: " + ratio);
	     }
		 ratio = temp_ratio.doubleValue();
		 return ratio;
	}

}
