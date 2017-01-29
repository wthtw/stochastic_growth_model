``` MATLAB
clear; clc;
tol = 0.0001; %Tolerance
maxit = 100; %max iterations
Qcalc = [0.2;0.2]; %initial guess 
Qsav = [0.2 0.2]; %will save values of Q for each iteration
[q1,q2] = meshgrid(0:0.1:1,0:0.1:1);
Q = [q1;q2];
SumQ = q1 + q2;

c = [0.6; 0.8]; eta = 1.6; e = -1/eta;
fval1 = SumQ.^e + e*(q1+q2).^(e-1) .* q1 - c(1)*q1;
fval2 = SumQ.^e + e*(q1+q2).^(e-1) .* q2 - c(2)*q2;

%Solve using the Newton Method
for it=1:maxit
    [fval,fjac] = cournot(Qcalc); %Use cournot.m function to determine fval and fjac
    Qcalc = Qcalc - fjac\fval;
    Qsav = [Qsav; Qcalc']; %save the values for Q found in that step
    if norm(fval) < tol, break, end
end

[C1,h1] = contour(q1,q2,fval1,[0 0]); %plot the countour where fval1 = 0
hold on
[C2,h2] = contour(q1,q2,fval2,[0 0]); %plot the countour where fval2 = 0
hold on
plot(Qsav(:,1),Qsav(:,2),':o'); %plot the points found through the iterations
colormap cool
```
