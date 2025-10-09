let S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

const BASES = ['A', 'G', 'C', 'T']

let DINUCLEOTIDES = []
let TRINUCLEOTIDES = []
let NUMBER_OF_DINUCLEOTIDES = {}
let NUMBER_OF_TRINUCLEOTIDES = {}
let TOTAL_NUM_DINUCLEOTIDES = S.length - 1
let TOTAL_NUM_TRINUCLEOTIDES = S.length - 2

for (let i = 0; i < BASES.length; i++) {
    for (let k = 0; k < BASES.length; k++) {
        DINUCLEOTIDES.push(`${BASES[i]}${BASES[k]}`)
    }
}

for (let i = 0; i < BASES.length; i++) {
    for (let j = 0; j < BASES.length; j++) {
        for (let k = 0; k < BASES.length; k++) {
            TRINUCLEOTIDES.push(`${BASES[i]}${BASES[j]}${BASES[k]}`)
        }
    }
}

console.log(DINUCLEOTIDES)
console.log(TRINUCLEOTIDES)

for (let d of DINUCLEOTIDES) {
    let count = 0
    for (let j = 0; j < S.length - 1; j++) {
        if (S.slice(j, j + 2) === d) count++
    }
    NUMBER_OF_DINUCLEOTIDES[d] = count
}

for (let t of TRINUCLEOTIDES) {
    let count = 0
    for (let j = 0; j < S.length - 2; j++) {
        if (S.slice(j, j + 3) === t) count++
    }
    NUMBER_OF_TRINUCLEOTIDES[t] = count
}

let totalDin = Object.values(NUMBER_OF_DINUCLEOTIDES).reduce((a, b) => a + b, 0)
let totalTri = Object.values(NUMBER_OF_TRINUCLEOTIDES).reduce((a, b) => a + b, 0)

console.log("Number of Dinucleotides:", totalDin)
console.log("Number of Trinucleotides:", totalTri)


for (let d of DINUCLEOTIDES) {
    let count = NUMBER_OF_DINUCLEOTIDES[d] || 0
    let pct = (count / TOTAL_NUM_DINUCLEOTIDES) * 100
    console.log(`${d}: ${pct}%`)
}

for (let t of TRINUCLEOTIDES) {
    let count = NUMBER_OF_TRINUCLEOTIDES[t] || 0
    let pct = (count / TOTAL_NUM_TRINUCLEOTIDES) * 100
    console.log(`${t}: ${pct}%`)
}
