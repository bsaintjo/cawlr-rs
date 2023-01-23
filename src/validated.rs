//! Data that has been validated for passing to training models

pub struct ValidSampleData(Vec<f64>);

impl ValidSampleData {
    pub fn validated(xs: Vec<f64>) -> Option<Self> {
        let xs: Vec<f64> = xs
            .into_iter()
            .filter(|&x| (x > 40.) && (x < 170.))
            .collect();
        if xs.len() < 2 {
            None
        } else {
            Some(ValidSampleData(xs))
        }
    }

    pub fn inner(self) -> Vec<f64> {
        self.0
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_validated() {
        let case = Vec::new();
        let xs = ValidSampleData::validated(case);
        assert!(xs.is_none(), "empty");

        let case = vec![f64::NAN; 10];
        let xs = ValidSampleData::validated(case);
        assert!(xs.is_none(), "NAN");

        let case = vec![f64::INFINITY; 10];
        let xs = ValidSampleData::validated(case);
        assert!(xs.is_none(), "INFINITY");

        let case = vec![f64::NEG_INFINITY; 100];
        let xs = ValidSampleData::validated(case);
        assert!(xs.is_none(), "NEG_INFINITY");

        let case = vec![100.0, 100.0, 140.0, 0.0, -0.0];
        let xs = ValidSampleData::validated(case);
        assert!(xs.is_some(), "hundreds");

        let case = vec![1.5028554895297472e233, 0.0];
        let xs = ValidSampleData::validated(case);
        assert!(xs.is_none(), "large finite");
    }
}
