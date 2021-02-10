% Programa para análise das bolhas (ou qq sejam as partículas) na
% folha de mica
%
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

ParticleNo{NumOfFiles} = zeros;
Area{NumOfFiles} = zeros; Mean{NumOfFiles} = zeros;
Min{NumOfFiles} = zeros; Max{NumOfFiles} = zeros;
Perim{NumOfFiles} = zeros; Circ{NumOfFiles} = zeros;
Feret{NumOfFiles} = zeros; Median{NumOfFiles} = zeros;
PercentArea{NumOfFiles} = zeros; FeretX{NumOfFiles} = zeros;
FeretY{NumOfFiles} = zeros; FeretAngle{NumOfFiles} = zeros;
MinFeret{NumOfFiles} = zeros; AR{NumOfFiles} = zeros;
Round{NumOfFiles} = zeros; Solidity{NumOfFiles} = zeros;
% Nova variável indicando o no. de part.
% em cada imagem
NumOfParticles(NumOfFiles) = zeros; StdDev{NumOfFiles} = zeros;
Label{NumOfFiles} = zeros; Radius{NumOfFiles} = zeros;
AltNano{NumOfFiles} = zeros;

%% Loop de importação do arquivo
% Cada coluna do arquivo .csv gerado pelo FIJI é uma variável
% ellipsis (...) indica quebra de linha

startRow = 2;
endRow = inf;

for i = 1:NumOfFiles
    
    filename{i} = strcat(path, name{i});
    [ParticleNo{i},Label{i},Area{i},Mean{i},StdDev{i},Min{i},Max{i},Perim{i},Circ{i},...
        Feret{i},Median{i},PercentArea{i},FeretX{i},FeretY{i},...
        FeretAngle{i},MinFeret{i},AR{i},Round{i},...
        Solidity{i}] = importfileFIJI(filename{i}, startRow, endRow);
    NumOfParticles(i) = length(ParticleNo{i});
    
    Radius{i} = (1/2)*sqrt(Area{i}/pi);
    AltNano{i} = Max{i}*1000000000;
    
end

%% Conversão p/ Força Lateral -> V -> nN
%{
conv = 65.66; % Fator de conversão do JPK

LateralForceMax{NumOfFiles} = zeros;
LFMean{NumOfFiles} = zeros; LFsig{NumOfFiles} = zeros;

for i=1:NumOfFiles
    
    LateralForceMax{i} = conv*Max{i};
    LFMean{i} = conv*Mean{i};
    LFsig{i} = conv*StdDev{i};
    
end
%}

%% Atribuindo cada arquivo à uma Força normla aplicada
%{
Set1 = 1:1:8;
Set2 = 0;
Setpoint = cat(2, Set1, Set2);
%}
Setpoint = 0:1:8;
k = 0.5; % N/m
Str = 56.84; % V/um
Adh = 3.533; % V

[TotalLoad, Fad, Fk] = LoadJPK(Str, Adh, k, Setpoint);
ApplLoad = Fk - Fad;
roundLoad = round(TotalLoad);

%% Ajuste linear p/ Altura-área
%{
p{NumOfFiles} = zeros;
S(NumOfFiles) = struct('R', [], 'df', [], 'norm', []);
Rsqr(NumOfFiles) = zeros; Radj(NumOfFiles) = zeros;
y{NumOfFiles} = zeros;
%
for i = 1:NumOfFiles
    
    [p{i}, S(i)] = polyfit(Area{i}, Max{i}, 1);
    %xs{i} = linspace(min(Area{i}, max(Area{i},200);
    y{i} = polyval(p{i}, Area{i});
    Rsqr(i) = 1 - S(i).norm^2 / norm(Max{i}-mean(Max{i}))^2;
    Radj(i) = 1 - (1 - Rsqr(i))*(NumOfParticles(i) - 1/...
        (NumOfParticles(i) - 2));
    
end
%}

pfit = struct;
S = struct;

for i = 1:NumOfFiles
    
    pfit(i).Load = roundLoad(i);
    [pfit(i).vec, pfit(i).S] = polyfit(Area{i}, AltNano{i}, 1);
    pfit(i).y = polyval(pfit(i).vec, Area{i});
    pfit(i).Rsqr = 1 - pfit(i).S.normr.^2 / norm(AltNano{i}-mean(AltNano{i})).^2;
    pfit(i).Radj = 1 - (1 - pfit(i).Rsqr)*(NumOfParticles(i) - 1)/...
        (NumOfParticles(i) - 2);
    
    
end

%% Correlações Área-Altura e Perímetro-Altura das partículas

close all

% Criação das legendas
legstr{NumOfFiles} = zeros;
figArea(NumOfFiles) = zeros; figPer(NumOfFiles) = zeros;
roundLoad(NumOfFiles) = zeros;

% Diretório onde os gráficos serão gravados
plotpath = strcat(path, 'Plots\');

% Loop para nomear arquivo
newstr{NumOfFiles} = zeros;
plotstr1{NumOfFiles} = zeros; plotstr2{NumOfFiles} = zeros;
legArea{NumOfFiles} = zeros;

for ii = 1:NumOfFiles
    
    newstr{ii} = name{ii}(1:end-4);
    
    plotstr1{ii} = strcat(plotpath, newstr{ii}, '-Area_Height_Corr.png');
    plotstr2{ii} = strcat(plotpath, newstr{ii}, '-Perim_Height_Corr.png');
    
end

% Plots
for k = 1:NumOfFiles
    
    figure(), figArea(k) = plot(Area{k}, AltNano{k}, 'bo');
    xlabel('Área (\mum^2)', 'FontSize', 16);
    ylabel('Altura (nm)', 'FontSize', 16);
    roundLoad(k) = round(TotalLoad(k));
    legstr{k} = strcat('F_N = ', {' '}, num2str(roundLoad(k)), {' '}, ' nN');
    legstr2 = 'Ajuste linear';
    
    hold on
    
    plot(Area{k}, pfit(k).y, 'r-', 'LineWidth', 2)
    legArea{k} = legend(legstr{k}, legstr2);
    set(legArea{k}, 'location', 'northwest');
    set(gca, 'fontsize', 14)
    saveas(gcf, plotstr1{k}, 'png');
    
    hold off
    
    figure(), figPer(k) = plot(Perim{k}, AltNano{k}, 'ks',...
        'MarkerFaceColor', 'k', 'MarkerSize', 8);
    grid on;
    xlabel('Perímetro (\mum)', 'FontSize', 16);
    ylabel('Altura (nm)', 'FontSize', 16);
    legend(legstr{k}, 'Location', 'NorthEast');
    set(gca, 'fontsize', 14)
    saveas(gcf, plotstr2{k}, 'png');
    
end

%% Histogramas

close all

nbins = 10;
%{
nelmPer{NumOfFiles} = zeros; nelmArea{NumOfFiles} = zeros;
nelmCirc{NumOfFiles} = zeros; nelmMax{NumOfFiles} = zeros;
centPer{NumOfFiles} = zeros; centArea{NumOfFiles} = zeros;
centCirc{NumOfFiles} = zeros; centMax{NumOfFiles} = zeros;
%}
histPer(NumOfFiles) = zeros; histArea(NumOfFiles) = zeros;
histCirc(NumOfFiles) = zeros; histAltNano(NumOfFiles) = zeros;
histAR(NumOfFiles) = zeros; histLF(NumOfFiles) = zeros;
histRadius(NumOfFiles) = zeros;
%}

% Diretório onde os gráficos serão gravados
histpath = strcat(path, 'Hists\');

% Loop para nomear arquivo
histstr1{NumOfFiles} = zeros; histstr2{NumOfFiles} = zeros;
histstr3{NumOfFiles} = zeros; histstr4{NumOfFiles} = zeros;
histstr5{NumOfFiles} = zeros; histstr6{NumOfFiles} = zeros;
histstr7{NumOfFiles} = zeros;

for m = 1:NumOfFiles
    
    histstr1{m} = strcat(histpath, newstr{m}, '-Perim-Hist.png');
    histstr2{m} = strcat(histpath, newstr{m}, '-Area-Hist.png');
    histstr3{m} = strcat(histpath, newstr{m}, '-Circ-Hist.png');
    histstr4{m} = strcat(histpath, newstr{m}, '-Max-Height-Hist.png');
    histstr5{m} = strcat(histpath, newstr{m}, '-Aspect-Ratio-Hist.png');
    histstr6{m} = strcat(histpath, newstr{m}, '-LF-Mean-Hist.png');
    histstr7{m} = strcat(histpath, newstr{m}, '-Radius-Hist.png');
end

% Histogramas
for j = 1:NumOfFiles
    
    [nelmPer, centPer] = hist(Perim{j}, nbins);
    [nelmArea, centArea] = hist(Area{j}, nbins);
    [nelmCirc, centCirc] = hist(Circ{j}, nbins);
    [nelmMax, centMax] = hist(AltNano{j}, nbins);
    [nelmAR, centAR] = hist(AR{j}, nbins);
    [nelmRad, centRad] = hist(Radius{j}, nbins);
    
    figure(), histPer(j) = bar(centPer, nelmPer);
    xlabel('Perímetro (\mum)', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    set(gca, 'fontsize', 14)
    saveas(gcf, histstr1{j}, 'png');
    
    figure(), histArea(j) = bar(centArea, nelmArea);
    xlabel('Área (\mum^2)', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    set(gca, 'fontsize', 14)
    saveas(gcf, histstr2{j}, 'png');
    
    figure(), histCirc(j) = bar(centCirc, nelmCirc);
    xlabel('Circularidade', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthWest');
    set(gca, 'fontsize', 14)
    saveas(gcf, histstr3{j}, 'png');
    
    figure(), histAltNano(j) = bar(centMax, nelmMax);
    xlabel('Altura (nm)', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    set(gca, 'fontsize', 14)
    saveas(gcf, histstr4{j}, 'png');
    
    figure(), histAR(j) = bar(centAR, nelmAR);
    xlabel('Razão de aspecto', 'FontSize', 16);
    ylabel('Frequência', 'FontSize', 16);
    legend(legstr{j}, 'Location', 'NorthEast');
    set(gca, 'fontsize', 14)
    saveas(gcf, histstr5{j}, 'png');
    
    figure(), histRadius(j) = bar(centRad, nelmRad);
    xlabel('Raio do círculo (\mum)', 'FontSize', 16)
    ylabel('Frequência', 'FontSize', 16)
    legend(legstr{j}, 'Location', 'NorthEast')
    set(gca, 'fontsize', 14)
    
end

%% Arquivo de resumo dos dados com media e desvio das variaveis

savename = strcat(path, 'Resumo Results w Erode.txt');

% Media das variaveis

muMax(NumOfFiles) = zeros; sigMax(NumOfFiles) = zeros;
muArea(NumOfFiles) = zeros; sigArea(NumOfFiles) = zeros;
muPer(NumOfFiles) = zeros; sigPer(NumOfFiles) = zeros;
muCirc(NumOfFiles) = zeros; sigCirc(NumOfFiles) = zeros;
muAR(NumOfFiles) = zeros; sigAR(NumOfFiles) = zeros;
muRad(NumOfFiles) = zeros; sigRad(NumOfFiles) = zeros;
%muLF(NumOfFiles) = zeros; sigLF(NumOfFiles) = zeros;

for jj = 1:NumOfFiles
    
    muMax(jj) = mean(Max{jj});
    sigMax(jj) = std(Max{jj});
    muArea(jj) = mean(Area{jj});
    sigArea(jj) = std(Area{jj});
    muPer(jj) = mean(Perim{jj});
    sigPer(jj) = std(Perim{jj});
    muCirc(jj) = mean(Circ{jj});
    sigCirc(jj) = std(Circ{jj});
    muAR(jj) = mean(AR{jj});
    sigAR(jj) = std(AR{jj});
    muRad(jj) = mean(Radius{jj});
    sigRad(jj) = std(Radius{jj});
    %muLF(jj) = mean(LFMean{jj});
    %sigLF(jj) = std(LFMean{jj});
    
end

% Para aquivos .jpk, em que a altura é dada em m, e não em nm
muAlt = 1000000000*muMax;
sigAlt = 1000000000*sigMax;

fid = fopen(savename, 'w');

% Criando cabeçalho do arquivo

h1 = 'Arquivo\tFn (nN)\tQuant. Part.\tAltura med. (nm)\tdesvio\t';
h2 = 'Area med. (um2)\tdesvio\tPerim. med. (um)\tdesvio\t';
h3 = 'Circ. med\tdesvio\tAspect Ratio\tdesvio\r\n';
header = strcat(h1,h2,h3);
fprintf(fid, header);
fclose(fid);

% Preenchendo arquivo até a penúltima linha
% Ordem a ser seguida:
% Arquivo; Força normal; Quant. de Part.; Altra media; desv. alt.; Area
% med; desv area; per med; desv per; circ med; desv circ; AR med; desv AR
fID = fopen(savename, 'a');

for kk=1:NumOfFiles-1
    
    fprintf(fID, '%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', ...
        newstr{kk},roundLoad(kk),NumOfParticles(kk),muAlt(kk),sigAlt(kk), ...
        muArea(kk),sigArea(kk),muPer(kk),sigPer(kk),muCirc(kk),sigCirc(kk),...
        muAR(kk),sigAR(kk));
    
end

fclose(fID);

% Preenche a última linha sem introduzir uma nova linha (\r\n)
% no arquivo .txt

FID = fopen(savename, 'a');

for kk=NumOfFiles:NumOfFiles
    
    fprintf(fID, '%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f', ...
        newstr{kk},roundLoad(kk),NumOfParticles(kk),muAlt(kk),sigAlt(kk), ...
        muArea(kk),sigArea(kk),muPer(kk),sigPer(kk),muCirc(kk),sigCirc(kk),...
        muAR(kk),sigAR(kk));
    
end

fclose(FID);

%% Cria um arquivo de resumo em csv
% Ainda em andamento...
%{
csvStr = strcat(path, 'Resumo Results.csv');

nrows = NumOfFiles;
ncols = 13;

A(nrows, ncols) = zeros;

for kk=1:ncols
    
    A(:,kk) = [newstr{kk},roundLoad(kk),NumOfParticles(kk),...
        muMax(kk),sigMax(kk), muArea(kk),sigArea(kk),muPer(kk),...
        sigPer(kk),muCirc(kk),sigCirc(kk),...
       muAR(kk),sigAR(kk)];
    
end

newheader = {'Arquivo', 'Fn (nN)','Quant. Part.', 'Altura med. (m)',...
    'desvio', 'Area med. (um2)', 'desvio', 'Perim. med. (um)',...
    'desvio', 'Circ. med', 'desvio', ...
    'Aspect Ratio', 'desvio'};
T = array2table(A);
T.Properties.Variablenames(1:end) = newheader;
writetable(T, csvStr);
%}
%% Plots das medias em fn da força normal

close all;

% Strings c/ nomes dos arqs.
strPer =  strcat(plotpath, 'Perimeter x load.png');
strArea =  strcat(plotpath, 'Area x load.png');
strMax =  strcat(plotpath, 'MaxHeight x load.png');
strCirc =  strcat(plotpath, 'Circ x load.png');
strAR =  strcat(plotpath, 'AR x load.png');
strNUM = strcat(plotpath, 'NumOfParticles x load.png');
strRadj = strcat(plotpath, 'Radj x load.png');
strRadius = strcat(plotpath, 'Radius x load.png');


% Plots
figure(), plot(roundLoad, muPer, 'ko',...
    'MarkerFaceColor', 'k', 'MarkerSize', 8);
ylabel('Perímetro (\mum)', 'FontSize', 16);
xlabel('Força normal (nN)', 'FontSize', 16);
set(gca, 'fontsize', 14)
grid on
saveas(gcf, strPer, 'png');

figure(), plot(roundLoad, muArea, 'ko',...
    'MarkerFaceColor', 'k', 'MarkerSize', 8);
ylabel('Área (\mum^2)', 'FontSize', 16);
xlabel('Força normal (nN)', 'FontSize', 16);
set(gca, 'fontsize', 14)
grid on
saveas(gcf, strArea, 'png');

figure(), plot(roundLoad, muMax, 'ko',...
    'MarkerFaceColor', 'k', 'MarkerSize', 8);
ylabel('Altura máxima (m)', 'FontSize', 16);
xlabel('Força normal (nN)', 'FontSize', 16);
set(gca, 'fontsize', 14)
grid on
saveas(gcf, strMax, 'png');

figure(), plot(roundLoad, muCirc, 'ko',...
    'MarkerFaceColor', 'k', 'MarkerSize', 8);
ylabel('Circularidade', 'FontSize', 16);
xlabel('Força normal (nN)', 'FontSize', 16);
set(gca, 'fontsize', 14)
grid on
saveas(gcf, strCirc, 'png');

figure(),plot(roundLoad, muAR, 'ko',...
    'MarkerFaceColor', 'k', 'MarkerSize', 8);
ylabel('Razão de aspecto', 'FontSize', 16);
xlabel('Força normal (nN)', 'FontSize', 16);
set(gca, 'fontsize', 14)
grid on
saveas(gcf, strAR, 'png');

figure(), plot(roundLoad, NumOfParticles, 'ko',...
    'MarkerFaceColor', 'k', 'MarkerSize', 8);
ylabel('No. de partículas segmentadas', 'FontSize', 16);
xlabel('Força normal (nN)', 'FontSize', 16);
ylim([57 67]);
set(gca, 'fontsize', 14)
grid on
saveas(gcf, strNUM, 'png');

figure() 
for i=1:NumOfFiles
    
    plot(pfit(i).Load, pfit(i).Radj, 'ks',...
    'MarkerFaceColor', 'k');
    hold on

end

ylim([0.6 0.7]);
xlim([25 105]);
xlabel('Força normal (nN)', 'FontSize', 16)
ylabel('R^2 ajustado', 'FontSize', 16)
set(gca, 'fontsize', 14)
saveas(gcf, strRadj, 'png');

hold off

figure(), errorbar(roundLoad, muRad, sigRad, 'ks',...
    'MarkerFaceColor', 'k');
xlabel('Força normal (nN)', 'FontSize', 16)
ylabel('Raio do círculo (\mum)', 'FontSize', 16)
set(gca, 'fontsize', 14)
saveas(gcf, strRadius, 'png')

%% Subplots de histogramas

close all;

% Subplot c/ Circularidade e Aspect Ratio

output_str2{NumOfFiles} = zeros;

for i=1:NumOfFiles
    
    output_str2{i} = strcat(histpath, 'Subplot\', ...
        'Perim-Area-Altura-', num2str(i), '.png');
    
    figure(), subplot(2,1,1)
    hist(Circ{i}, 10);
    xlabel('Circularidade', 'FontSize', 14);
    ylabel('Frequência', 'FontSize', 14);
    ylim([0 30]);
    set(gca, 'ytick', 0:10:30)
    set(gca, 'fontsize', 12)
    legend(legstr{i}, 'Location', 'NorthWest')
    
    subplot(2,1,2)
    hist(AR{i}, 10);
    ylim([0 35])
    xlabel('Razão de aspecto', 'FontSize', 14);
    ylabel('Frequência', 'FontSize', 14);
    legend(legstr{i});
    set(gca, 'fontsize', 12)
    saveas(gcf, output_str2{i}, 'png')
    
end

% Subplot c/ Area, Perim e Altura

output_str1{NumOfFiles} = zeros;

for i=1:NumOfFiles
    
    output_str1{i} = strcat(histpath, 'Subplot\', ...
        'Perim-Area-Altura-', num2str(i), '.png');
    
    figure(), subplot(2,2,1)
    hist(Perim{i}, 10);
    ylim([0 20])
    xlabel('Perímetro (\mum)', 'FontSize', 14);
    ylabel('Frequência', 'FontSize', 14);
    legend(legstr{i}, 'Location', 'NorthEast');
    set(gca, 'fontsize', 10)
    
    subplot(2,2,2)
    hist(Area{i}, 10);
    ylim([0 20])
    xlabel('Área (\mum^2)', 'FontSize', 14);
    ylabel('Frequência', 'FontSize', 14);
    legend(legstr{i}, 'Location', 'NorthEast')
    
    subplot(2,2,3:4)
    hist(AltNano{i}, 10);
    ylim([0 15])
    xlabel('Altura (nm)', 'FontSize', 14);
    ylabel('Frequência', 'FontSize', 14);
    legend(legstr{i});
    set(gca, 'fontsize', 12)
    saveas(gcf, output_str1{i}, 'png')
    
end

%% Raio de curvatura
%{
curv{NumOfFiles} = zeros; % Raio de curvatura
microAlt{NumOfFiles} = zeros; % Altura em micrometros
muCurv(NumOfFiles) = zeros; stdCurv(NumOfFiles) = zeros;
muMicroAlt(NumOfFiles) = zeros;

for i=1:NumOfFiles
    
    microAlt{i} = 1000000*Max{i};
    muMicroAlt(i) = mean(microAlt{i});
    curv{i} = sqrt(2*Radius{i}.*microAlt{i} - microAlt{i}.^2);
    muCurv(i) = mean(curv{i});
    stdCurv(i) = std(curv{i});
end

muCurvAll = mean(muCurv(:));
stdCurvAll = std(muCurv(:));
muMicroAltAll = mean(muMicroAlt(:));
stdMicroAltAll = std(muMicroAlt(:));
%}
%% histogramas raio de curvatura
%{
for i=1:NumOfFiles
    
    figure(), hist(curv{i}, nbins);
    xlabel('Raio de curvatura (\mum)', 'FontSize', 16)
    ylabel('Frequência', 'FontSize', 16)
    set(gca, 'fontsize', 12)
    legend(legstr{i}, 'Location', 'NorthEast')
    
end
%}
