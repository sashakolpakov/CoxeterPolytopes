FROM ubuntu:20.04

ARG ALVIN_COMMIT=d92006a842cb4668e7c89e8c8e2d30337e4b0f26

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        bison \
        build-essential \
        ca-certificates \
        cmake \
        flex \
        git \
        graphviz \
        libgmp-dev \
        libpcre3-dev \
        libpari-dev \
        libtool \
        libxml2-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/rgugliel/AlVin.git /opt/AlVin \
    && cd /opt/AlVin \
    && git checkout --detach "${ALVIN_COMMIT}" \
    && git submodule update --init --recursive

COPY vendor.sha256 /tmp/alvin-vendor.sha256
RUN cd /opt/AlVin/lib \
    && sha256sum -c /tmp/alvin-vendor.sha256 \
    && tar -xzf crlibm-1.0beta5.tgz \
    && tar -xzf gaol-4.2.0.tar.gz \
    && tar -xzf igraph-0.7.1.tar.gz \
    && tar -xjf eigen-eigen-323c052e1731.tar.bz2

COPY gaol-aarch64.patch /tmp/gaol-aarch64.patch
RUN cd /opt/AlVin/lib/gaol-4.2.0 \
    && patch -p1 < /tmp/gaol-aarch64.patch \
    && autoreconf -fi

WORKDIR /opt/AlVin/lib/crlibm
RUN ./configure && make -j1 && make install && ldconfig

WORKDIR /opt/AlVin/lib/gaol-4.2.0
RUN ./configure --disable-simd --disable-preserve-rounding --disable-asm \
        --disable-verbose-mode --with-mathlib=crlibm \
    && make -j1 \
    && make install \
    && ldconfig

WORKDIR /opt/AlVin/lib/igraph-0.7.1
RUN ./configure && make -j1 && make install && ldconfig

WORKDIR /opt/AlVin/lib/eigen-eigen-323c052e1731
RUN mkdir build && cd build && cmake .. && make install

WORKDIR /opt/AlVin/build
COPY alvin-norm-powers.patch /tmp/alvin-norm-powers.patch
COPY alvin-aarch64.patch /tmp/alvin-aarch64.patch
COPY alvin-rc7-exact-bounds.patch /tmp/alvin-rc7-exact-bounds.patch
COPY alvin-rc7-angle14.patch /tmp/alvin-rc7-angle14.patch
COPY alvin-rc7-fraction-conjugation.patch /tmp/alvin-rc7-fraction-conjugation.patch
COPY alvin-rc7-spatial-r1.patch /tmp/alvin-rc7-spatial-r1.patch
RUN cd /opt/AlVin && patch -p1 < /tmp/alvin-norm-powers.patch
RUN cd /opt/AlVin && patch -p1 < /tmp/alvin-aarch64.patch
RUN cd /opt/AlVin && patch -p1 < /tmp/alvin-rc7-exact-bounds.patch
RUN cd /opt/AlVin && patch -p1 < /tmp/alvin-rc7-angle14.patch
RUN cd /opt/AlVin && patch -p1 < /tmp/alvin-rc7-fraction-conjugation.patch
RUN cd /opt/AlVin && patch -p1 < /tmp/alvin-rc7-spatial-r1.patch
RUN cmake .. && make -j1

ENTRYPOINT ["/opt/AlVin/build/alvin"]
