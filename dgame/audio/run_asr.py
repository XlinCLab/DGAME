import argparse
import whisper


def transcribe_audio(file_path: str, model_name: str, language: str):
    # Load model
    model = whisper.load_model(model_name)

    # Run transcription
    result = model.transcribe(
        file_path,
        language=language,
        task="transcribe",
        # # Avoids aggressive context smoothing that may reduce disfluencies
        # condition_on_previous_text=False,
    )

    return result["text"]


def main():
    parser = argparse.ArgumentParser(description="Perform ASR transcription with OpenAI Whisper.")

    parser.add_argument(
        "audio_path",
        type=str,
        help="Path to the input .wav file"
    )

    parser.add_argument(
        "--model",
        type=str,
        default="medium",
        help="Whisper model size (tiny, base, small, medium, large)"
    )

    parser.add_argument(
        "--language",
        type=str,
        default="de",
        help="Language code (default: de for German)"
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional path to save transcript as .txt file"
    )

    args = parser.parse_args()

    text = transcribe_audio(
        file_path=args.audio_path,
        model_name=args.model,
        language=args.language
    )

    print("\n--- TRANSCRIPTION ---\n")
    print(text)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"\nSaved transcript to: {args.output}")


if __name__ == "__main__":
    main()