from typing import Tuple

import pynq


class FpgaContext:
    # Translates memory name to attribute name of overlay.
    _MEM_TABLE = {ord(c): None for c in '[]'}

    def __init__(self, bitstream: str):
        self.shard_count = 0
        self.interval_count = 0
        self._overlay = pynq.Overlay(bitstream)
        self.edges = []
        self.vertices = []
        for arg in self._overlay.SSSP.args.values():
            if arg.name == f'edges_{self.shard_count}':
                mem = getattr(self._overlay, arg.mem.translate(self._MEM_TABLE))
                self.edges.append(
                    pynq.allocate(
                        shape=250 * 1024 * 1024 // 8,
                        dtype='i4, f4',
                        target=mem,
                    ))
                self.shard_count += 1
            elif arg.name == f'vertices_{self.interval_count}':
                mem = getattr(self._overlay, arg.mem.translate(self._MEM_TABLE))
                self.vertices.append(
                    pynq.allocate(
                        shape=(250 * 1024 * 1024 // 16,),
                        dtype='i4, f4, i4, i4',
                        target=mem,
                    ))
                self.interval_count += 1
        self._metadata = pynq.allocate(
            shape=(128 * 128 // 8,),
            dtype='i8',
            target=self._overlay.PLRAM0,
        )
        self._heap_array = pynq.allocate(
            shape=(8388608,),
            dtype='u8, u8, u8, u8',
            target=self._overlay.DDR0,
        )
        self._heap_index = pynq.allocate(
            shape=(256 * 1024 * 1024 // 8,),
            dtype='u8',
            target=self._overlay.HBM17,
        )

    def call(
        self,
        vertex_count: int,
        root: int,
        root_offset: int,
        root_degree: int,
    ) -> None:
        """Call the FPGA kernel.

        Edges must be sync_to_device() before calling this function.
        """
        for buf in self.vertices:
            buf[:vertex_count // self.interval_count].sync_to_device()

        self._overlay.SSSP.call(
            vertex_count,  # vertex_count
            root,
            root,
            root,
            root_offset,
            root_degree,
            self._metadata,
            *self.edges,
            *self.vertices,
            self._heap_array,
            self._heap_index,
        )

        for buf in self.vertices:
            buf[:vertex_count // self.interval_count].sync_from_device()


fpga_context = FpgaContext(
    '/home/blaok/tmp/tapa-sssp/1cdd1bb0cb/build/SSSP.xilinx_u280_xdma_201920_3.hw.xclbin'
)
