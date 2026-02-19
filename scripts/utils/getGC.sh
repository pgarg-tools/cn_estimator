awk ' \
BEGIN { \
    FS=""; \
    cg=0; \
    t=0; \
} \
{ \
    if ($1 != ">") { \
        for (i = 1; i <= NF; i++) { \
            if ($i ~ /[ACTGactg]/) { \
                t++;
            } \
            if ($i ~ /[CGcg]/) { \
                cg++;
            } \
        } \
    } \
    else { \
        if (t > 0) { \
            gc_perc = 1000*(cg/t); \
            if (gc_perc % 10 < 5){ gc_perc = int(gc_perc/10); } else {gc_perc = 1 + int(gc_perc/10); } \
	        printf("%s\t%d\t%d\t%f\t%.0f\n",h,cg,t,cg/t,gc_perc); \
            cg = 0; \
            t = 0; \
            gc_perc = "NA"; \
        } \
        h = substr($0,2); \
    } \
} \
END { \
    gc_perc = 1000*(cg/t); \
    if(gc_perc % 10 < 5){ gc_perc = int(gc_perc/10); } else {gc_perc = 1 + int(gc_perc/10); } \
	printf("%s\t%d\t%d\t%f\t%.0f\n",h,cg,t,cg/t,gc_perc); \
}' $1
