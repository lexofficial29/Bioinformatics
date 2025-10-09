let S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

let DINUCLEOTIDES = []
let TRINUCLEOTIDES = []

for (let i = 0; i < S.length - 1; i++) {
    DINUCLEOTIDES.push(`${S.at(i)}${S.at(i+1)}`)
}

for (let i = 0; i < S.length - 2; i++) {
    TRINUCLEOTIDES.push(`${S.at(i)}${S.at(i+1)}${S.at(i+2)}`)
}

console.log(DINUCLEOTIDES.length)
console.log(TRINUCLEOTIDES.length)