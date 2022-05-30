objs=""

echo "CC = g++"
echo "LD = g++"
echo "CC_FLAGS = -std=c++11 -O2 -I. -I./include/ -I/usr/include/freetype2/"
echo "LD_FLAGS = "
echo "LIBS = -lglut -lGLEW -lglfw3 -lGLU -lGL -lX11 -ldl -lXi -lXrandr -lXxf86vm -lXinerama -lXcursor -lrt -lm -pthread -lfreetype -lAntTweakBar"
echo ""
echo ".PHONY: all"
echo "all: cs467"
echo ""

for arg in "$@"
do
  sh generate.sh $arg
  name=`echo $arg | sed "s/\///g" | tr [a-z] [A-Z]`
  objs="$objs \$($name)"
done

echo "OBJS =$objs\n"

dirs=""

for arg in "$@"
do
  name=`echo $arg | sed "s/\///g"`
  echo "${name}_dir:"
  dirs="$dirs ${name}_dir"
  echo "\tmkdir -p out/$name\n"
done

echo ".PHONY: directories"

echo "directories: $dirs\n"

echo "cs467: directories \$(OBJS)"
echo "\t\$(LD) \$(LD_FLAGS) \$(OBJS) \$(LIBS) -o out/cs467"
echo ""

echo "clean:"
echo "\trm -r out/*"
