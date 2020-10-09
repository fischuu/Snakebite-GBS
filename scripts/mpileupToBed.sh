awk 'BEGIN{
    for(r=getline; r>0;){
        for(s=e=$2; (r=getline)>0 && $2<=e+1; e=$2);
        print s==e ? s : $1"\t"s"\t"e
    }
    exit -r
}' $1
