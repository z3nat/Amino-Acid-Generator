#!/usr/bin/perl
#amino.pl
use strict;
use warnings;


my %bondValenceCharges = (1 => 1, 2 => (-2), 3 => (-3), 4 => (-4));
my %atomMwt = (
1 => 1.0079,
2 => 15.9994,
3 => 14.0067,
4 => 12.0107
); 

my %formula_MWt = (
 89.0929 =>'C3H7NO2 Alanine',
 121.1579 =>'C3H7NO2S Cystenine',
 133.1024 => 'C4H7NO4 Aspartic Acid',
 147.1289 => 'C5H9NO4 Glutamic Acid',
 165.1887 => 'C9H11NO2 Phenylalanine',
 75.0664 => 'C2H5NO2 Glycine',
 155.1542 =>'C6H9N3O2 Histidine',
 131.1724 => 'C6H13NO2 Isoleucine',
 146.1870 => 'C6H14N2O2 Lysine',
 131.1724 => 'C6H13NO2 Leucine',
 149.2109 => 'C5H11NO2S Methionine',
 132.1176 => 'C4H8N2O3 Asparagine',
 115.1301 => 'C5H9NO2 Proline',
 146.1441 => 'C5H10N2O3 Glutamine',
 174.2004 => 'C6H14N4O2 Arginine',
 105.0923 => 'C3H7NO3 Serine',
 119.1188 => 'C4H9NO3 Threonine',
 117.1459 => 'C5H11N02 Valine',
 204.2247 => 'C11H12N2O2 Tryptophan',
 181.1881 => 'C9H11NO3 Tyrosine')
;
my @AA_MWt = (89.0929,121.1579,133.1024,147.1289,165.1887,75.0664,155.1542,131.1724,146.1870,131.1724,149.2109,132.1176,115.1301,146.1441,174.2004,105.0923,119.1188,117.1459,204.2247,181.1881);
my @similar_MWt =(
);
#chose to use all 20 amino acid weights instead of 17 or 18, because in the future I could adjust my program to account for isomers.
my $range = 4;
my $minimum = 1;
my $CountOfRuns = 0;
my $z=0;
my $valence = 0 ;
my $molecular_weight=0; 
my $AA = 0; 
my $Hcount =0;
my $Ocount =0 ;
my $Ncount = 0;
my $Ccount =0; 
my $diff;
my $i;
my $l =0;
my $j=0;
my @difference;
my @sort_difference= sort {$a <=> $b} @difference;
my $closest_match= 0;

until($AA != 0 && $molecular_weight <= 204.2247 && $valence == 0){ 
	
	do{ $z = int(rand($range)) + $minimum;				
				
		if ($valence < 0 and $z > 1){
			
			$valence+=2;
			} #when bond valence is - and random number valence is -
			
		if ($z == 1){
			
			$valence +=$bondValenceCharges{1};
			
			$molecular_weight += $atomMwt{1};
			
			$Hcount++;
			
		}elsif ($z == 2){
			
			$valence +=$bondValenceCharges{2};
			
			$molecular_weight += $atomMwt{2};
			
			$Ocount++;
			
		}elsif ($z == 3){
			
			$valence += $bondValenceCharges{3};
			
			$molecular_weight += $atomMwt{3};
			
			$Ncount++;
			
		}elsif ($z== 4){
			
			$valence += $bondValenceCharges {4};
			
			$molecular_weight += $atomMwt{4};
			
			$Ccount++;
			
		}if ($Ccount != 0 && $z == 4 && $Ccount%3== 0){ #adding double bond here for every 3rd carbon
			
			$valence+= 2;
			
		}if ($valence < -4){ #adding 3 hydrogens when valence falls below -4
			
			$valence+=3;		
			
		}if (exists($formula_MWt{$molecular_weight}) == 1 && $valence == 0){ #adding an amino acid to AA count
			
			$AA++;
			
			print "**Amino Acid Created!**\n";
			
			print "Amino Acid is: $formula_MWt{$molecular_weight}\n";
		}
	
		 $CountOfRuns++;
		
		if ($valence == 0 && $molecular_weight< 204.2247 && $AA != 0){
			
			push(@similar_MWt,$molecular_weight);	
			
			for($l .. $#AA_MWt){
				
				$diff = $AA_MWt[$l] - $molecular_weight; #finding the difference b/w each stable compound and amino acid
		
				$l++;
					
					if($diff>0){ #if the difference is positive, adding it to an array
						
						push(@difference, $diff)
					}
			}				
			@sort_difference = sort {$a <=> $b} @difference;
		
			for ($j.. $#similar_MWt){
			
				$closest_match = $similar_MWt[$j] + $sort_difference[0];
			
			 	$j++;			
			}
			if(exists($formula_MWt{$closest_match}) == 1){
				
				print "Closest match was to $formula_MWt{$closest_match} with a molecular weight difference of $sort_difference[0]\n";
			}
			
		}
		
		if ($CountOfRuns%10 == 0 ||$valence == 0 || $AA ==1 ||exists($formula_MWt{$closest_match})==1 ){ #printing here for every 10th run or when a stable compound is formed or when a amino acid is formed
			
			 print "Current Number of Open Bonds is: $valence\n"; 
			
			 print "Current Compound weight is : $molecular_weight \n";
			
			 print "Current Compound fomula is: C$Ccount H$Hcount N$Ncount O$Ocount \n";
			
			 print "Amino Acid Count is: $AA\n";
			
			 print "Count of random trials: $CountOfRuns\n";
		
			 #sleep 1;	
		}
		
		if($CountOfRuns%10 == 0 && exists($formula_MWt{$closest_match})!=1 && $AA == 0 ){
			
			print"There is no current closest match \n";
			
			print "\n";
		}
		if( exists($formula_MWt{$closest_match})!=1 && $valence == 0){
			
			print "There is no current closest match\n";
			
			print "\n";
		} 
	}
	while ($valence != 0 && $molecular_weight <= 204.2247 && $AA == 0); 
		
		if($molecular_weight >= 204.2247){
			
			$valence = 0;
			
			$Hcount = 0;
			
			$Ocount = 0;
			
			$Ncount = 0;
			
			$Ccount = 0;
			
			$molecular_weight = 0;
		}
	
}
print "*******************************************************\n";

print "Program done!\n" ; 
exit;
#Student Assignment Submission Form 
#==================================
#I declare that the attached assignment is wholly my/our
#own work in accordance with Seneca Academic Policy.  No part of this 
#assignment has been copied manually or electronically from any 
#other source (including web sites) or distributed to other students. 
#Name: Zena Teferi Student ID#: 109682179
#-------------------------------------------------------------
