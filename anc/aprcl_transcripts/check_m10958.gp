p = readstr("prp_m10958.txt")[1];
n = eval(p);
print("Loaded n with ", #digits(n), " digits.");
gettime();
result = isprime(n, 2);
ms = gettime();
write("aprcl_done_m10958.txt", result);
print("RESULT m=10958: ", result, " (gettime=", ms, "ms)");
quit;
