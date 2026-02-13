
__device__ __forceinline__ char iupac_decode(u8 iupac) {
  switch (iupac) {
    case 0b0001: return 'A';
    case 0b0010: return 'C';
    case 0b0100: return 'G';
    case 0b1000: return 'T';
    case 0b0101: return 'R';
    case 0b1010: return 'Y';
    case 0b0110: return 'S';
    case 0b1001: return 'W';
    case 0b1100: return 'K';
    case 0b0011: return 'M';
    case 0b1110: return 'B';
    case 0b1101: return 'D';
    case 0b1011: return 'H';
    case 0b0111: return 'V';
    case 0b1111: return 'N';
  }
  printf("invalid iupac encoding!");
  return '?';
}

__device__ __forceinline__ bool iupac_match(u8 a, u8 b) {
  return (a & b) != 0;
}
