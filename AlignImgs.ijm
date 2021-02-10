// Alinhamento dos canais ida e volta de topografia

run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/MF-18-MultiLayer-2019.12.10-10.50.48-FWD.txt]");
run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/MF-18-MultiLayer-2019.12.10-10.50.48-BWD.txt]");
selectWindow("MF-18-MultiLayer-2019.12.10-10.50.48-FWD.txt");
run("Images to Stack", "name=Stack title=[] use keep");
run("StackReg", "transformation=Translation");


// Alinhamento manual das imagens é dado pelo comando Translate
// Image > Transform > Translate

run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Test Matrix 2.txt]");
run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Test Matrix 1.txt]");
selectWindow("Test Matrix 1.txt");
run("Translate...", "x=9 y=0 interpolation=None");
imageCalculator("Subtract create 32-bit", "Test Matrix 2.txt","Test Matrix 1.txt");
selectWindow("Result of Test Matrix 2.txt");