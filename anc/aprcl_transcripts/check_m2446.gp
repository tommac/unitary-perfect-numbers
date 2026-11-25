p = readstr("prp_m2446.txt")[1];
n = eval(p);
print("Loaded n with ", #digits(n), " digits.");
gettime();
result = isprime(n, 2);
ms = gettime();
write("aprcl_done_m2446.txt", result);
print("RESULT m=2446 (368-digit): ", result, " (gettime=", ms, "ms)");
quit;
