
syms x1 x2 x3 real;

%% Συνάρτηση.
%f = x1^3 - 3*x1^2 + x2^2;
%f = x1^3+x1*x3^2+3*x1^2+x2^2+2*x3^2;
%f = x1^3+x2^2+x3^2-3*x1;
%f = x1^4+x1^2-6*x1*x2+3*x2^2;
%f = x1*x2^2 + x1^3*x2 - x1*x2;
%f = x1^2+2*x2^2+5*x3^2-2*x1*x2-4*x2*x3-2*x3;
f=x1^4+2*x1^2*x2^2+x2^4;

%f = -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) - exp( 0.5*(cos(2*pi*x1)+cos(2*pi*x2) ) )+ exp(1) +20;
%f =  (1.5 -x1 +x1*x2)^2 + (2.25 - x1 + x1*x2^2)^2 + (2.625 - x1 + x1*x2^3)^2;
%f = (x1^2 + 2*x2^2 + 3*x3^2)*exp(-(x1^2 + x2^2 + x3^2));
%f = [1 + (x1 + x2 +1)^2 * (19 -14*x1 + 3*x1^2 -14*x2 + 6*x1*x2 + 3*x2^2)] * [30 + (2*x1 - 3*x2)^2 * (18 -32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)];
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

D_values = zeros(length(x),1); %Αρικοποίηση του πίνακα  για τις υπο-ελάσσονες τιμές με μηδενικά.
for i=1:rows %Για όσα είναι τα κρίσιμα σημεία.
    HF_star = subs(HF, x, X(i,:)); %%Αντικατέστησε τις τιμές του κρίσιμου σημείου στον Εσσιανό πίνακα.
    for j=1:length(x) %Για όσες είναι οι ανεξάρτητες μεταβλητές.
        D = HF_star(1:j,1:j); %Φτιάξε τον πίνακα κάθε φορά για μια ελλάσονα ορίζουσα.
        D_values(j)= det(D); %Υπολόγισε την ορίζουσα του πίνακα και καταχωρησέ την στον πίνακα με τις τιμές των οριζουσών.
    end
    
    %Περιπτώσεις.
    total_flag = 1; %Δήλωση ότι δεν έχει εκτυπωθεί κάποιο μήνυμα.
    %1)
    flag = 1; %Έστω ότι όλες οι τιμές είναι θετικές.
    for k=1:length(D_values)
        if D_values(k) <= 0  %Αν κάποια δεν είναι τότε 
            flag = 0;   %δεν είναι τοπικός ελαχιστοποιητής.
            break;
        end
    end
    if flag && total_flag
        fprintf("[");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local minimizer.\n");
        total_flag = 0;
    end
    
    %2)
    flag = 1; %Έστω ότι οι τιμες σε όλες τις περιττές θέσεις είναι αρνητικές και σε όλες τις άρτιες θετικές.
    for k=1:length(D_values)
        if (mod(k,2)==1 && D_values(k)>0) || (mod(k,2)==0 && D_values(k)<0 ) %Αν κάποιο στοιχείο σε περιττή θέση
            flag = 0; %είναι θετικό ή κάποιο στοιχείο σε άρτια θέση είναι αρνητικό τότε δεν είναι τοπικός μεγιστοποιητής.
            break;
        end
    end
    if flag && total_flag
        fprintf("[");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local maximizer.\n");
        total_flag = 0;
    end
    
    %3)
    flag = 1; %Έστω όλα τα στοιχεία σε άρτιες θέσεις είναι θετικά.
    for k=1:length(D_values)
        if mod(k,2)==0 && D_values(k) < 0 %Αν κάποιο στοιχείο σε άρτια θέση είναι αρνητικό 
            flag = 0;  %τότε δεν είναι τοπικός βελτιστοποιητής.
            break;
        end
    end
    if flag == 0 && total_flag
        fprintf("[");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is not local optimizer.\n");
        total_flag = 0;
    end
    
    %4)
    flag_zero = 0; %Έστω ότι δεν υπάρχουν μηδενικές τιμές.
    flag_neg = 0; %Έστω ότι δεν υπάρχουν αρνητικές τιμές.
    flag_pos = 0; %Έστω ότι δεν υπάρχουν θετικές τιμές.
    
    %Έλεγχος στοιχείων (Το πρώτο γίνεται ξεχωριστά). Αν βρεθεί 0 ή θετικό ή
    %αρνητικό στοιχείο, η αντίστοιχη μεταβλητή γίνεται ίση με την μονάδα.
    prev = D_values(1);
    if  prev > 0
        flag_pos = 1;
    elseif prev < 0
        flag_neg = 1;
    else
        flag_zero = 1;
    end
    if flag_zero == 0
        for k=2:length(D_values)
            if D_values(k) == 0
                flag_zero = 1;
                break;
            end
            if D_values(k)/prev > 0
                flag_pos = 1;
            else
                flag_neg = 1;
            end
            prev = D_values(k);
        end
    end
    if (flag_zero == 0 && flag_pos && flag_neg && total_flag) %Αν δεν υπάρχουν μηδενικά και υπάρχουν θετικά και
        fprintf("[");                                         %αρνητικά στοιχεία, τότε δεν είναι τοπικός βελτιστοποιητής.
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is not local optimizer case.\n");
        total_flag = 0;
    end
    
    if total_flag %Άμα δεν συμβαίνει τίποτα από τα προηγούμενα τότε ο αλγόριθμος δεν μπορεί να βγάλει κάποιο συμπέρασμα.
        fprintf("Can't figure if [");
        for k=1:columns-1
            fprintf("%s,",X(i,k));
        end
        fprintf("%s",X(i,columns));
        fprintf("]");
        fprintf(" is local optimizer.\n");
    end

end