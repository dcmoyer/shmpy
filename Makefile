

default: third_party shmpy

spherical_harmonics/bin/spherical_harmonics.a:
	$(MAKE) spherical_harmonics -C spherical_harmonics

third_party: spherical_harmonics/bin/spherical_harmonics.a

shmpy: third_party
	python shmpy/setup.py build_ext --build-lib shmpy/ --force
	#python shmpy/setup.py build_ext --force

