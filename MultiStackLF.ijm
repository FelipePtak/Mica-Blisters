run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/MF-18-MultiLayer-2019.12.10-10.50.48-FWD.txt]");
run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Topography/MF-18-MultiLayer-2019.12.10-10.50.48-BWD.txt]");
selectWindow("MF-18-MultiLayer-2019.12.10-10.50.48-FWD.txt");
run("Images to Stack", "name=Stack-Topo title=[] use");
run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Lateral Force/WSxM/MF-18-MultiLayer-2019.12.10-10.50.48-FWD.txt]");
run("Text Image... ", "open=[F:/Documentos/Felipe/PUC/DOUTORADO/Mica/MF-18/Setpoint 2/Text Matrices/Lateral Force/WSxM/MF-18-MultiLayer-2019.12.10-10.50.48-BWD.txt]");
selectWindow("MF-18-MultiLayer-2019.12.10-10.50.48-FWD.txt");
run("Images to Stack", "name=Stack-LF title=[] use");
selectWindow("Stack-Topo");
run("StackReg", "transformation=Translation");
run("MultiStackReg", "stack_1=Stack-Topo action_1=[Use as Reference] file_1=[] stack_2=Stack-LF action_2=[Align to First Stack] file_2=[] transformation=Translation");

// Salvando a matriz de transf.

run("MultiStackReg", "stack_1=Stack-Topo action_1=[Use as Reference] file_1=[] stack_2=Stack-LF action_2=[Align to First Stack] file_2=[] transformation=Translation save");