if [ "$#" -ne 1 ]; then
  echo "usage: ./generate.sh <directory>"
  exit 1
fi

dirname=`echo $1 | sed "s/\///g" | tr [a-z] [A-Z]`
objs=""
build_rules=""

for i in `find $1 -name "*.cpp"`
do
  o_file=out/`echo $i | sed "s/.cpp/.o/g"`
  objs="$objs $o_file"
  rule="$o_file: $i"
  build="\$(CC) \$(CC_FLAGS) -c -o $o_file $i\n"
  build_rules="$build_rules\n$rule\n\t$build"
done

echo "$dirname =$objs"
echo $build_rules
