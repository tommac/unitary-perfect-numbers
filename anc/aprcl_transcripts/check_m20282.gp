p = readstr("prp_m20282.txt")[1];
n = eval(p);
print("Loaded n with ", #digits(n), " digits.");
gettime();
result = isprime(n, 2);
ms = gettime();
write("aprcl_done_m20282.txt", result);
print("RESULT m=20282: ", result, " (gettime=", ms, "ms)");
quit;
