p = readstr("prp_m17398.txt")[1];
n = eval(p);
print("Loaded n with ", #digits(n), " digits.");
gettime();
result = isprime(n, 2);
ms = gettime();
write("aprcl_done_m17398.txt", result);
print("RESULT m=17398: ", result, " (gettime=", ms, "ms)");
quit;
