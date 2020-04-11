% Test gemaakt om resultaten te kunnen vergelijken met c++
VB = 5;VH = 6; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
Varray = 1/2*ones(VB,VH); %Volledig metaal materiaal
Q = 200;
BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
BC1 = [['D',1,1,293];['D',2,VH-1,293];['D',VH,VH,293]]; % Rechter 
BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
BC3 = [['D',1,1,273];['D',2,VH-1,273];['D',VH,VH,273]];% Linker 
[Sol,K] = Harmonic_FVM(VB,VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3,2);

 SOL = reshape(Sol,[VB,VH]);
 figure(1); surf(SOL); 
 title("Testcase");
