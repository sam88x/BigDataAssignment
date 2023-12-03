DATA bigData;
INFILE 'C:\Users\PhotonUser\My Files\OneDrive\Files\STAT 4110\fat.dat.txt';
INPUT Case 3-5 Brozek 10-13 Siri 18-21 Density 24-29 Age 36-37 Weight 40-45 Height 49-53 Adiposity 58-61 FatFree 65-69 Neck 74-77 Chest 81-85 Abdomen 89-93 Hip 97-101 Thigh 106-109 Knee 114-117 Ankle 122-125 Biceps 130-133 Forearm 138-141 Wrist 146-149;
RUN;

proc print data=bigData;
run;

proc export data=bigData
	outfile = 'C:\Users\PhotonUser\My Files\OneDrive\Files\STAT 4110\bigData.csv'
	dbms = csv
	replace;
run;
