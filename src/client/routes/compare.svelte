<script lang="ts">
	import state from '$meta';

	const { query } = state;

	const ordering = query.datasets.split(',');
</script>

{#await fetch(`/plots?datasets=${query.datasets}&gene=${query.gene}&groupBy=${query.groupBy}&splitBy=${query.splitBy}`).then((res) => res.json())}
	<i class="fa-solid fa-spinner"></i> Generating Plots...
{:then plots}
	<div class="row">
		{#each ordering as dataset}
			<div class="col">
				<h3 class="dataset">{dataset.replaceAll('_', ' ')}</h3>
				<img src={plots[dataset].clustering} alt="{dataset} clustering" />
				<img src={plots[dataset].violin} alt="{dataset} violin" />
			</div>
		{/each}
	</div>
{:catch err}
	<h1>Error: {err}</h1>
{/await}

<style lang="scss">
	@keyframes spin {
		0% {
			transform: rotate(0deg);
		}

		100% {
			transform: rotate(360deg);
		}
	}

	img {
		width: 300px;
	}

	.dataset {
		text-align: center;
	}

	:global(.fa-spinner) {
		animation: spin 1s linear infinite;
	}
</style>
