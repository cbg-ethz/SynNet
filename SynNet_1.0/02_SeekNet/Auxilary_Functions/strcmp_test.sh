function strcmp_test
{
if [ "$1" = "$2" ] ; then
return 0
fi
if expr "$1" ">" "$2" ; then 
return 1
fi
return -1
}
