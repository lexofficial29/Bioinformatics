const SEQUENCE = "ATTTCGCCGATA"

let alphabet = {}
let frequency = {}
let total = 0

for(let i = 0 ; i < SEQUENCE.length; i++) {
    if(isNaN(alphabet[SEQUENCE.at(i)])) {
        alphabet[SEQUENCE.at(i)] = 1
        total += 1
    }
    else {
        alphabet[SEQUENCE.at(i)] += 1
        total += 1
    }
}

for (let key in alphabet) {
    frequency[key] = ((alphabet[key] / total) * 100).toFixed(2) + "%"
}

console.log(frequency)


