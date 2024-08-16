# Geometric Algebra Library Written In Rust

```
	let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
	let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

	let out = SeaMultivector4D::<f32>::from_str("3,1,0>-e3 +e13 +e14 +e23 +e123 +e124").unwrap();
	assert_eq!(
		&mv_a * &mv_b,
		out
	);

	let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
	let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

	let out = SeaMultivector4D::<f32>::from_str("3,1,0>-e3 +e23").unwrap();
	assert_eq!(
		&mv_a | &mv_b,
		out
	);

	let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
	let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

	let out = SeaMultivector4D::<f32>::from_str("3,1,0>e13 +e14 +e123 +e124").unwrap();
	assert_eq!(
		&mv_a ^ &mv_b,
		out
	);

	let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
	let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

	let out = SeaMultivector4D::<f32>::from_str("3,1,0>-e3 +e23").unwrap();
	assert_eq!(
		&mv_a << &mv_b,
		out
	);

	let mv_a = SeaMultivector4D::<f32>::from_str("3,1,0>e1 + e12").unwrap();
	let mv_b = SeaMultivector4D::<f32>::from_str("3,1,0>e3 + e123 + e4").unwrap();

	let out = SeaMultivector4D::<f32>::from_str("3,1,0>").unwrap();
	assert_eq!(
		&mv_a >> &mv_b,
		out
	);
```