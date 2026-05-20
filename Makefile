RAPTOR_DIR=/home/wqzhang/opt/raptor
LLVM_VER=20

default: a.out

a.out: main.cpp
	$(RAPTOR_DIR)/bin/raptor-clang++ -stdlib=libc++ -O3 -flto=full -fuse-ld=lld \
            -Rpass=raptor \
            -o $@ $<  \
            -Xlinker "--load-pass-plugin=$(RAPTOR_DIR)/lib/LLDRaptor-$(LLVM_VER).so" \
            -L$(RAPTOR_DIR)/lib -lRaptor-RT-$(LLVM_VER) -lmpfr -lstdc++

#a.out: main.cpp
#	$(RAPTOR_DIR)/bin/raptor-clang++ -stdlib=libc++ -O3 -Rpass=raptor -o $@ $<

clean:
	rm -f *.o

distclean: clean
	rm -f a.out
