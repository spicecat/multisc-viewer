<script lang="ts">
import type { Snippet } from "svelte";

const {
	tags = $bindable(),
	tag,
	children,
}: { tags?: string[]; tag?: string; children: Snippet } = $props();
</script>

<button
	type="button"
	class="chip preset-tonal"
	onclick={() => {
		if (tags && tag && !tags.includes(tag)) tags.push(tag);
	}}>{@render children()}</button
>
