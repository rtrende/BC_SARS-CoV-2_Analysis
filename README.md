**Protocol to Analyze BC SARS-CoV-2 NGS Results**

*CURRENTLY SOME DETAILS ARE SPECIFIC TO BOON LAB, PROTOCOL TO BE UPDATED FOR MORE GENERAL ANALYSIS*

**Downloading Data and Setting Up Folders**

1. Make a new empty folder with your experiment name
2. Within that folder:
    1. Make empty folders called “barcode_counts” and “fastq_files”
    2. Copy the “ref” folder from a previous run (this just contains the reference SARS-CoV-2 genome)
    3. Copy the BC_analysis_BCs.py, BC_analysis_WT.py, and BC_analysis_part2.py files from a previous run
    4. Cut the bbmap folder from a previous run and paste into new folder
        1. Large folder, copying takes much longer and takes up more space on your hard drive
    5. By the end, your folder should look like this:

![Example Folder Structure](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture1.png?raw=true)

1. Open the windows command line (cmd)
2. Enter “bash” to run Ubuntu
3. Change directories into your experiment’s directory, then into fastq_files
4. Download the data with the following command:
    1. wget -r --no-parent \[copy URL from Jess and ML\]
    2. You should get a bunch of output that looks like this over and over:

![Example Script Output](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture3.png?raw=true)

1. The data will be buried several folders deep (e.g. for my example it was in fastq_files/ dsildata.htcf.wustl.edu/Boon-10997548/200-768). Cut all files out of this deep folder and move them to fastq_files
2. Delete dsildata.htcf.wustl.edu folder

**Counting the BCs in each sample**

1. Change directories into your experiment’s directory
    1. If you’re still in fastq_files, cd ../
2. Run the following series of commands (part2.py should run relatively fast, but the other two steps could take 1-3 hours depending on sequencing depth)
    1. python3 BC_analysis_BCs.py
    2. python3 BC_analysis_part2.py
        1. After running this, open giant_panda_clean.csv in Excel and re-save it as \[your_experiment_name\]\_raw.xlsx
    3. python3 BC_analysis_WT.py
    4. python3 BC_analysis_part2.py
        1. There is a chance you get an error on this last one (if any samples do not have anything map to the WT sequence, it makes part2 angry, and most experiments have at least 1 sample with no WT sequence)… if so:
            1. Open \[your_experiment_name\]\_raw.xlsx in Excel
            2. Insert a new row under the headers
            3. Call row “WT_\_D614G” (🡨 Two underscores, the name needs to be 9 letters long to match the length of the BC names to make Matlab happy)
            4. Navigate to the barcode_counts folder
            5. For each sample, look at the file that ends in .barcode_counts.txt and see how many reads have “TTTTCGCTT” (making the preview pane visible in file explorer makes this much faster)
            6. Add that information to the WT_\_D614G row of your spreadsheet
        2. Example:
            1. Add a new row with WT_\_D614G

![Adding_D614G_Step1](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture4.png?raw=true)

- - - 1. Look at the .barcode_counts.txt folder for the first sample, see that 11 reads map to TTTTCGCTT

![Adding_D614G_Step2](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture5.png?raw=true)
![Adding_D614G_Step3](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture6.png?raw=true)


- - - 1. Rinse and repeat for all subsequent samples. For samples with no reads, enter 0; for samples that have reads mapping to something other than TTTTCGCTT, ignore all other entries (example below)

![Adding_D614G_Step4](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture7.png?raw=true)
![Adding_D614G_Step5](https://github.com/rtrende/BC_SARS-CoV-2_Analysis/blob/main/images/Picture8.png?raw=true)

**Checking column names**

1. Check that all columns have an experiment number (containing T#, where # can be any number)
2. Check that all tissues collected are properly specified in each excel column. The list of tissues my current Matlab script can automatically ID is shown in the table below (case sensitive).
    1. Can add more at any time, just give me a heads up

| **Tissue** | **Abbreviation Matlab looks for** |     | **Tissue** | **Abbreviation Matlab looks for** |
| --- | --- | --- | --- | --- |
| Nasal turbinate | NT  |     | Right caudal lobe | RL3 |
| Trachea | Trach |     | Infracardial lobe | RL4 |
| Whole lung | WL  |     | Heart | Heart |
| Left lobe | LL  |     |     |     |
| Right apical lobe | RL1 |     |     |     |
| Right middle lobe | RL2 |     |     |     |

1. Check that all information you want matlab to know about is in the name of each excel column. For example, if you collect the donors at different times, make sure that info is in the file name. The list of things my current Matlab script can automatically detect, see the table below
    1. Note that you don’t have to include everything from this table below, just the variables you want Matlab to account for. For example, don’t need duration of exposure in the column name if that’s not something you varied

| Variable | What phrase does Matlab look for? | Example of acceptable column name |
| --- | --- | --- |
| Duration of exposure | hr  | T46_H09_8hr_NT |
| Route of exposure | AB, DC, Fom | T41_H11_AB_Trach |
| Time of collection | If donor, hpi<br><br>If contact, hpe | T41_H01_32hpi_WL<br><br>T43_H23_48hpe_LL |
| Animal sex | M or F | T42_H03_M_RL1 |

\*: donors and contacts are detected by animal number, where if animal number is less than or equal to the number of donors, that animal is marked as a donor. Otherwise it is marked as a contact

1. Check that all the information you want Matlab to be able to detect is surrounded by underscores
2. Copy the excel file with the WT_\_D614G data and properly labeled columns into a new data analysis folder for Matlab
3. Open Matlab
4. Run Trim_BC_Data and follow the prompts (you should just have to input the file name. Be sure to include the file extension!)
