const SEQUENCE = "ATTTCGCCGATA"

let alphabet = {}

for(let i = 0 ; i < SEQUENCE.length; i++) {
    if(isNaN(alphabet[SEQUENCE.at(i)])) {
        alphabet[SEQUENCE.at(i)] = 1
    }
    else {
        alphabet[SEQUENCE.at(i)] += 1
    }
}

console.log(alphabet)