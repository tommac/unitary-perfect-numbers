p = readstr("prp_m10294.txt")[1];
n = eval(p);
print("Loaded n with ", #digits(n), " digits.");
gettime();
result = isprime(n, 2);
ms = gettime();
write("aprcl_done_m10294.txt", result);
print("RESULT m=10294: ", result, " (gettime=", ms, "ms)");
quit;
