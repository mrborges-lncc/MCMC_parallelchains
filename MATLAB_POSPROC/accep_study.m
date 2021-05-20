clear;
close all;

dados = load('../twoStage/data.dat');
hist(dados(:,2),40);