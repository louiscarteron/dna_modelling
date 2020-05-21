class Error {

  constructor(stages) {
    this.errors = stages.reduce((acc, val) => {
      acc[val] = {
        'inserts': 0,
        'deletes': 0,
        'substitutions': 0
      };
      return acc;
    }, {});

    this.incrementError = (stage, errorType, value = 1) => {
      this.errors[stage][errorType] += value;
    }

    this.printReport = () => {
      console.log(this.errors);
    }

  }


}

export default Error;