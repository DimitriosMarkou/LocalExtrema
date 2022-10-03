
syms x1 x2 x3 real;

%% Συνάρτηση.
%f = x1^3 - 3*x1^2 + x2^2;
%f = x1^3+x1*x3^2+3*x1^2+x2^2+2*x3^2;
%f = x1^3+x2^2+x3^2-3*x1;
%f = x1^4+x1^2-6*x1*x2+3*x2^2;
%f = x1^2 - 6*x1*x2 + 2*x2^2 + 10*x1 + 2*x2 -5;
%f = (x1^2 + 2*x2^2 + 3*x3^2)*exp(-(x1^2 + x2^2 + x3^2));
%f = (x1^2 + 2*x2^2 + 3*x3^2)*exp(-(x1^2 + x2^2 + x3^2));

%f =  (1.5 -x1 +x1*x2)^2 + (2.25 - x1 + x1*x2^2)^2 + (2.625 - x1 + x1*x2^3)^2;
%f = sin(x1+x2) + (x1-x2)^2 -1.5*x1 +2.5*x2 + 1;
x = [x1, x2];

Gradf = sym(zeros(length(x),1)); %Αρχικοποίηση πίνακα με μηδενικά.
for i=1:length(x)
    Gradf(i,1) = diff(f,x(1,i)); %Έυρεση πρώτης παραγώγου για κάθε μεταβλητή.
end

HF = hessian(f,x) %Υπολογισμός του Εσσιανού πίνακα.

X = table2array(struct2table(solve(Gradf==0,x))); %Εύρεση κρίσιμων σημείων.

[rows, columns] = size(X);

for i=1:rows
    for j=1:columns
        X(i,j)=vpa(X(i,j)); %Αναγκάζει να υπολογιστεί η συνολική πράξη. 
    end
end

for i=1:rows %Για όσες είναι οι γραμμές του Χ (το σύνολο των κρίσιμων σημείων).
    HF_star = subs(HF, x, X(i,:)); %Αντικατέστησε τις τιμές του κρίσιμου σημείου στον Εσσιανό πίνακα.
    eig_star=eig(HF_star); %Υπολόγισε την ιδιοτιμή του Εσσιανού πίνακα.
  
    
    %Περιπτώσεις.
    %1)
    flag_pos = 0; %Έστω ότι δεν υπάρχουν μηδενικές τιμές.
    flag_neg = 0; %Έστω ότι δεν υπάρχουν αρνητικές τιμές.
    flag_zero = 0; %Έστω ότι δεν υπάρχουν μηδενικές τιμές.
    %Αν βρεθεί 0 ή θετικό ή αρνητικό στοιχείο, η αντίστοιχη μεταβλητή γίνεται ίση με την μονάδα.
    for j=1:length(eig_star)
        eig_star = simplify(eig_star);
        if  cast(eig_star(j),'double') > 0
            flag_pos = 1;
        elseif cast(eig_star(j),'double') < 0
            flag_neg = 1;
        else 
            flag_zero = 1;
        end
    end
    
    if flag_pos && flag_neg==0 && flag_zero==0 %Αν είναι ολα θετικά τότε το σημείο έιναι τοπικό ελάχιστο.
        fprintf("[");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local minimizer.\n");
        
    elseif flag_neg && flag_pos==0 && flag_zero==0 %Αν είναι ολα αρνητικά τότε το σημείο έιναι τοπικό μέγιστο.
        fprintf("[");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local maximizer.\n");
  
    elseif flag_pos && flag_neg %Αν υπάρχουν και θετικά και αρνητικά στοιχεία τότε το σημείο δεν είναι τοπικός βελτιστοποιητής.
        fprintf("[");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is not local optimizer.\n");
        
    else %Για οποιαδήποτε άλλη περίπτωση ο αλγόριθμος δεν μπορεί να βγάλει κάποιο συμπέρασμα.
        fprintf("Can't figure if [");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local optimizer.\n");
    end
end