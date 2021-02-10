run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/MF-18-MultiLayer-2019.12.10-11.11.18-FWD.txt]");
run("Duplicate...", " ");
setOption("ScaleConversions", true);
run("8-bit");
run("Auto Threshold", "method=Otsu white");
run("Watershed");
run("Set Scale...", "distance=256 known=5 unit=µm global");
run("Set Measurements...", "area mean min perimeter shape feret's median area_fraction redirect=MF-18-MultiLayer-2019.12.10-11.11.18-FWD.txt decimal=3");
run("Analyze Particles...", "size=0.02-Infinity show=[Overlay Masks] display exclude clear");
saveAs("PNG", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Segmented Images/MF-18-MultiLayer-2019.12.10-11.11.18-FWD-1.png");
run("Distribution...", "parameter=Max or=21 and=0-0");
run("Distribution...", "parameter=Perim. or=21 and=0-0");
run("Distribution...", "parameter=Area or=21 and=0-0");
run("Distribution...", "parameter=Circ. or=21 and=0-0");
run("Distribution...", "parameter=AR or=21 and=0-0");

// Salvando resutados em .txt

run("Input/Output...", "jpeg=85 gif=-1 file=.txt use_file copy_row save_column save_row");

saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Results/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Results Otsu Seg.txt");
//selectWindow("Perim. Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Per Distr FIJI.txt");
//selectWindow("Max Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Max Distr FIJI.txt");
//selectWindow("Area Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Area Distr FIJI.txt");
//selectWindow("Circ. Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Circ Distr FIJI.txt");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD AspRat Distr FIJI.txt");


// Salvando resutados em .csv

run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row");

saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Results/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Results Otsu Seg.csv");
//selectWindow("Perim. Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Per Distr FIJI.csv");
//selectWindow("Max Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Max Distr FIJI.csv");
//selectWindow("Area Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Area Distr FIJI.csv");
//selectWindow("Circ. Distribution");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD Circ Distr FIJI.csv");
//saveAs("Results", "F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/Data Analysis/Distributions/MF-18-MultiLayer-2019.12.10-11.11.18-FWD AspRat Distr FIJI.csv");

run("Clear Results");
close();
close();
close();
close();
close();
close();
close();