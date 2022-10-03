syms x1 x2 x3 real;

%% Συναρτήσεις.
%f = x1^3 - 3*x1^2 + x2^2;
%f= x1^3 + x1*x3^2 + 3*x1^2 + x2^2 + 2*x3^2;
%f = x1^3 + x2^2 + x3^2 - 3*x1;
%f = x1^4 + x1^2 - 6*x1*x2 + 3*x2^2;
%f = 3*x1^4 + 3*x1^2*x2 - x2^3;
%f = -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) - exp( 0.5*(cos(2*pi*x1)+cos(2*pi*x2) ) )+ exp(1) +20;
%f =  (1.5 -x1 +x1*x2)^2 + (2.25 - x1 + x1*x2^2)^2 + (2.625 - x1 + x1*x2^3)^2;
f = 2*x1^2 - 1.05*x1^4 + x1^6/6 + x1*x2 + x2^2;
%f = sin(x1+x2) + (x1-x2)^2 -1.5*x1 +2.5*x2 + 1;
x = [x1, x2];

Gradf = sym(zeros(length(x),1)); %Αρχικοποίηση πίνακα με μηδενικά.
for i=1:length(x)
    Gradf(i,1) = diff(f,x(1,i)); %Έυρεση πρώτης παραγώγου για κάθε μεταβλητή.
end

HF = hessian(f,x); %Υπολογισμός του Εσσιανού πίνακα.


X = table2array(struct2table(solve(Gradf==0,x))); %Εύρεση κρίσιμων σημείων. (μετατοπή απο struct σε table)
[rows, columns] = size(X);
for i=1:rows
    for j=1:columns
        X(i,j)=vpa(X(i,j)); %Αναγκάζει να υπολογιστεί η συνολική πράξη.
    end
end

for i=1:rows   %Για όσες είναι οι γραμμές του Χ (το σύνολο των κρίσιμων σημείων).
    HF_star = subs(HF, x, X(i,:)); %Αντικατέστησε τις τιμές του κρίσιμου σημείου στον Εσσιανό πίνακα.
    Qf_star = x * HF_star * x.'; %Υπολόγισε το quadratic form του κρίσιμπου σημείου.
    Qf_star = simplify(Qf_star); %Aπλοποιήση του quadratic form.
 
    assume(x > 0 | x < 0)
    if isAlways(Qf_star > 0,'Unknown','false') %Αν είναι πάντα θετικό είναι τοπικός ελαχιστοποιητής.
        fprintf("Qf_star(");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf(")(");
        for k=1:columns-1
            fprintf("%s,",x(1,k));
        end
        fprintf("%s",x(1,columns));
        fprintf(") = %s > 0 so [",Qf_star);
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local minimizer.\n");
    
    elseif isAlways(Qf_star < 0,'Unknown','false') %Αν είναι πάντα αρνητικό είναι τοπικός μεγιστοποιητής.
        fprintf("Qf_star(");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf(")(");
        for k=1:columns-1
            fprintf("%s,",x(1,k));
        end
        fprintf("%s",x(1,columns));
        fprintf(") = %s < 0 so [",Qf_star);
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local maximizer.\n");
    
        %Αν δεν είναι ούτε μεγαλύτερο ή ίσο του μηδενός και ούτε μικρότερο ή ίσο του μηδενός τότε δεν είναι τοπικός βελτιστοποιητής.
    elseif isAlways(Qf_star >= 0,'Unknown','false') == 0 && isAlways(Qf_star <= 0,'Unknown','false') == 0 
        fprintf("Qf_star(");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf(")(");
        for k=1:columns-1
            fprintf("%s,",x(1,k));
        end
        fprintf("%s",x(1,columns));
        fprintf(") = %s > 0 and < 0 so [",Qf_star);
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is not local optimizer.\n");
    else
        fprintf("Qf_star(");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf(")(");
        for k=1:columns-1
            fprintf("%s,",x(1,k));
        end
        fprintf("%s",x(1,columns));
        fprintf(") = %s is semi-defined so can't figure if [",Qf_star);
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is  local optimizer.\n");
    end
    
end


