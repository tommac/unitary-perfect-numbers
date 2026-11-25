p = readstr("prp_m19066.txt")[1];
n = eval(p);
print("Loaded n with ", #digits(n), " digits.");
gettime();
result = isprime(n, 2);
ms = gettime();
write("aprcl_done_m19066.txt", result);
print("RESULT m=19066: ", result, " (gettime=", ms, "ms)");
quit;
