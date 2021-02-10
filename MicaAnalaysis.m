% Programa para análise das bolhas (ou qq sejam as partículas) na folha de
% mica
% Amostra: MF-18; Conjunto: Setpoint 2;
% Diretório: /DOUTORADO/Mica/MF-18/Setpoint 2/
% Arquivos: MF-18-Multilayer-2019.12.10-time.jpk

clear;
clc;
close all;

%% Pegando nome dos arquivos selecionando-os por GUI

[name, path] = uigetfile('*.csv', 'Multiselect', 'on');

%% Inicializando variaveis p/ otimização

NumOfFiles = length(name);
filename{NumOfFiles} = zeros;

ParticleNumber{NumOfFiles} = zeros;
Area{NumOfFiles} = zeros; Mean{NumOfFiles} = zeros;
Min{NumOfFiles} = zeros; Max{NumOfFiles} = zeros;
Perim{NumOfFiles} = zeros; Circ{NumOfFiles} = zeros;
Feret{NumOfFiles} = zeros; Median{NumOfFiles} = zeros;
PercentArea{NumOfFiles} = zeros; FeretX{NumOfFiles} = zeros;
FeretY{NumOfFiles} = zeros; FeretAngle{NumOfFiles} = zeros;
MinFeret{NumOfFiles} = zeros; AR{NumOfFiles} = zeros;
Round{NumOfFiles} = zeros; Solidity{NumOfFiles} = zeros;
NumOfParticles(NumOfFiles) = zeros; % Nova variável indicando o no. de part.
% em cada imagem

%% Loop de importação do arquivo
% Cada coluna do arquivo .csv gerado pelo FIJI é uma variável
% ellipsis (...) indica quebra de linha

startRow = 2;
endRow = inf;

for i = 1:NumOfFiles
    
    filename{i} = strcat(path, name{i});
    [ParticleNumber{i},Area{i},Mean{i},Min{i},Max{i},Perim{i},Circ{i},...
        Feret{i},Median{i},PercentArea{i},FeretX{i},FeretY{i},...
        FeretAngle{i},MinFeret{i},AR{i},Round{i},...
    Solidity{i}] = importFIJIResults(filename{i}, startRow, endRow);
    NumOfParticles(i) = length(ParticleNumber{i});
    
end

%% Atribuindo cada arquivo à uma Força normla aplicada

Set1 = 1:1:8;
Set2 = 0;
Setpoint = cat(2, Set1, Set2);
k = 0.5; % N/m
S = 56.84; % V/um
Adh = 3.533; % V

[TotalLoad, Fad, Fk] = LoadJPK(S, Adh, k, Setpoint);
ApplLoad = Fk - Fad;

%% Correlações Área-Altura e Perímetro-Altura das partículas

% Criação das legendas
legstr{NumOfFiles} = zeros;
figArea(NumOfFiles) = zeros; figPer(NumOfFiles) = zeros;
roundLoad(NumOfFiles) = zeros;

% Diretório onde os gráficos serão gravados
plotpath = strcat(path, 'Plots\');

% Loop para nomear arquivo
newstr{NumOfFiles} = zeros;
plotstr1{NumOfFiles} = zeros; plotstr2{NumOfFiles} = zeros;

for ii = 1:NumOfFiles
    
    newstr{ii} = name{ii}(1:end-21);
    
    plotstr1{ii} = strcat(plotpath, newstr{ii}, '-Area_Height_Corr.png');
    plotstr2{ii} = strcat(plotpath, newstr{ii}, '-Perim_Height_Corr.png');
    
end

% Plots
for k = 1:NumOfFiles
    
    figure(), figArea(k) = plot(Area{k}, Max{k}, 'bo');
    xlabel('Área (\mum^2)', 'FontSize', 16);
    ylabel('Altura (m)', 'FontSize', 16);
    roundLoad(k) = round(TotalLoad(k));
    legstr{k} = strcat('F_N = ', {' '}, num2str(roundLoad(k)), {' '}, ' nN');
    legend(legstr{k}, 'Location', 'NorthWest');
    saveas(gcf, plotstr1{k}, 'png');
    
    figure(), figPer(k) = plot(Perim{k}, Max{k}, 'ks', 'MarkerFaceColor',...
        'k', 'MarkerSize', 8);
    grid on;
    xlabel('Perímetro (\mum)', 'FontSize', 16);
    ylabel('Altura (m)', 'FontSize', 16);
    legend(legstr{k}, 'Location', 'NorthWest');
    saveas(gcf, plotstr2{k}, 'png');
    
end

%% Histogramas

close all

nbins = 21;
%{
nelmPer{NumOfFiles} = zeros; nelmArea{NumOfFiles} = zeros;
nelmCirc{NumOfFiles} = zeros; nelmMax{NumOfFiles} = zeros;
centPer{NumOfFiles} = zeros; centArea{NumOfFiles} = zeros;
centCirc{NumOfFiles} = zeros; centMax{NumOfFiles} = zeros;
%}
histPer(NumOfFiles) = zeros; histArea(NumOfFiles) = zeros;
histCirc(NumOfFiles) = zeros; histMax(NumOfFiles) = zeros;
histAR(NumOfFiles) = zeros;
%}

% Diretório onde os gráficos serão gravados
histpath = strcat(path, 'Hists\');

% Loop para nomear arquivo
histstr1{NumOfFiles} = zeros; histstr2{NumOfFiles} = zeros;
histstr3{NumOfFiles} = zeros; histstr4{NumOfFiles} = zeros;
histstr5{NumOfFiles} = zeros;

for m = 1:NumOfFiles
    
    histstr1{m} = strcat(histpath, newstr{m}, '-Perim-Hist.png');
    histstr2{m} = strcat(histpath, newstr{m}, '-Area-Hist.png');
    histstr3{m} = strcat(histpath, newstr{m}, '-Circ-Hist.png');
    histstr4{m} = strcat(histpath, newstr{m}, '-Max-Height-Hist.png');
    histstr5{m} = strcat(histpath, newstr{m}, '-Aspect-Ratio-Hist.png');
end

% Histogramas
for j = 1:NumOfFiles
    
    [nelmPer, centPer] = hist(Perim{j}, nbins);
    [nelmArea, centArea] = hist(Area{j}, nbins);
    [nelmCirc, centCirc] = hist(Circ{j}, nbins);
    [nelmMax, centMax] = hist(Max{j}, nbins);
    [nelmAR, centAR] = hist(AR{j}, nbins);
    
    figure(), histPer(j) = bar(centPer, nelmPer);
    xlabel('Perímetro (m)', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    saveas(gcf, histstr1{j}, 'png');
    
    figure(), histArea(j) = bar(centArea, nelmArea);
    xlabel('Área (\mum^2)', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    saveas(gcf, histstr2{j}, 'png');
    
    figure(), histCirc(j) = bar(centCirc, nelmCirc);
    xlabel('Circularidade', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthWest');
    saveas(gcf, histstr2{j}, 'png');
    
    figure(), histMax(j) = bar(centMax, nelmMax);
    xlabel('Altura (m)', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    saveas(gcf, histstr4{j}, 'png');
    
    figure(), histAR(j) = bar(centAR, nelmAR);
    xlabel('Razão de aspecto', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'Best');
    saveas(gcf, histstr5{j}, 'png');
    
end

%% Arquivo de resumo dos dados com media e desvio das variaveis

savename = strcat(path, 'Resumo Results.txt');

% Media das variaveis

muMax(NumOfFiles) = zeros; sigMax(NumOfFiles) = zeros;
muArea(NumOfFiles) = zeros; sigArea(NumOfFiles) = zeros;
muPer(NumOfFiles) = zeros; sigPer(NumOfFiles) = zeros;
muCirc(NumOfFiles) = zeros; sigCirc(NumOfFiles) = zeros;
muAR(NumOfFiles) = zeros; sigAR(NumOfFiles) = zeros;

for jj = 1:NumOfFiles
    
    muMax(jj) = mean(Max{jj});
    sigMax(jj) = std(Max{jj});
    muArea(jj) = mean(Area{jj});
    sigArea(jj) = std(Area{jj});
    muPer(jj) = mean(Per{jj});
    sigPer(jj) = std(Per{jj});
    muCirc(jj) = mean(Circ{jj});
    sigCirc(jj) = std(Circ{jj});
    muAR(jj) = mean(AR{jj});
    sigAR(jj) = std(AR{jj});
    
end

fid = fopen(savename, 'w');

% Criando cabeçalho do arquivo

header = 'Arquivo\tFn (nN)\tQuant. Partículas\tAltura med. (m)\tdesvio\tArea med. (um2)\tdesvio\tPerim. med. (m)\tdesvio\tCirc. med\tdesvio\tAspect Ratio\r\n';
fprintf(fid, header);
fclose(fid);

% Preenchendo arquivo até a penúltima linha
% Ordem a ser seguida:
% Arquivo; Força normal; Quant. de Part.; Altra media; desv. alt.; Area
% med; desv area; per med; desv per; circ med; desv circ; AR med; desv AR
fID = fopen(savename, 'a');

for kk=1:NumOfFiles-1
    
   fprintf(fID, '%s\t%f\t%f\t%f\t%f\t%f\%f\%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', ...
       newstr{kk},roundLoad(kk),NumOfParticles(kk),muMax(kk),sigMax(kk), ...
       muArea(kk),sigArea(kk),muPer(kk),sigPer(kk),muCirc(kk),sigCirc(kk),...
       muAR(kk),sigAR(kk));
   
end

fclose(fID);

% Preenche a última linha sem introduzir uma nova linha (\r\n)
% no arquivo .txt

FID = fopen(savename, 'a');

for kk=NumOfFiles:NumOfFiles
    
   fprintf(fID, '%s\t%f\t%f\t%f\t%f\t%f\%f\%f\t%f\t%f\t%f\t%f\t%f\t%f', ...
       newstr{kk},roundLoad(kk),NumOfParticles(kk),muMax(kk),sigMax(kk), ...
       muArea(kk),sigArea(kk),muPer(kk),sigPer(kk),muCirc(kk),sigCirc(kk),...
       muAR(kk),sigAR(kk));
   
end

fclose(FID);